/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceExtendedCustomCVForce.h"

#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/ReferenceTabulatedFunction.h"
#include "lepton/CustomFunction.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include "lepton/Operation.h"

using namespace OpenMMLab;
using namespace OpenMM;
using namespace Lepton;
using namespace std;

// This class allows us to update tabulated functions without having to recompile expressions
// that use them.
class ReferenceExtendedCustomCVForce::TabulatedFunctionWrapper : public CustomFunction {
public:
    TabulatedFunctionWrapper(vector<Lepton::CustomFunction*>& tabulatedFunctions, int index) :
            tabulatedFunctions(tabulatedFunctions), index(index) {
    }
    int getNumArguments() const {
        return tabulatedFunctions[index]->getNumArguments();
    }
    double evaluate(const double* arguments) const {
        return tabulatedFunctions[index]->evaluate(arguments);
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        return tabulatedFunctions[index]->evaluateDerivative(arguments, derivOrder);
    }
    CustomFunction* clone() const {
        return new TabulatedFunctionWrapper(tabulatedFunctions, index);
    }
private:
    vector<Lepton::CustomFunction*>& tabulatedFunctions;
    int index;
};

ReferenceExtendedCustomCVForce::ReferenceExtendedCustomCVForce(const ExtendedCustomCVForce& force) {
    int numCVs = force.getNumCollectiveVariables();
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < numCVs; i++)
        variableNames.push_back(force.getCollectiveVariableName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        paramDerivNames.push_back(force.getEnergyParameterDerivativeName(i));

    // Create custom functions for the tabulated functions.

    map<string, CustomFunction*> functions;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions(), NULL);
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        tabulatedFunctions[i] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
        functions[force.getTabulatedFunctionName(i)] = new TabulatedFunctionWrapper(tabulatedFunctions, i);
    }

    // Create the expressions.

    ParsedExpression energyExpr = Parser::parse(force.getEnergyFunction(), functions).optimize();
    energyExpression = energyExpr.createCompiledExpression();
    variableDerivExpressions.clear();
    for (auto& name : variableNames)
        variableDerivExpressions.push_back(energyExpr.differentiate(name).createCompiledExpression());
    paramDerivExpressions.clear();
    for (auto& name : paramDerivNames)
        paramDerivExpressions.push_back(energyExpr.differentiate(name).createCompiledExpression());
    globalValues.resize(globalParameterNames.size());
    cvValues.resize(numCVs);
    map<string, double*> variableLocations;
    for (int i = 0; i < globalParameterNames.size(); i++)
        variableLocations[globalParameterNames[i]] = &globalValues[i];
    for (int i = 0; i < numCVs; i++)
        variableLocations[variableNames[i]] = &cvValues[i];
    energyExpression.setVariableLocations(variableLocations);
    for (CompiledExpression& expr : variableDerivExpressions)
        expr.setVariableLocations(variableLocations);
    for (CompiledExpression& expr : paramDerivExpressions)
        expr.setVariableLocations(variableLocations);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

void ReferenceExtendedCustomCVForce::updateTabulatedFunctions(const ExtendedCustomCVForce& force) {
    // Create custom functions for the tabulated functions.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        if (tabulatedFunctions[i] != NULL) {
            delete tabulatedFunctions[i];
            tabulatedFunctions[i] = NULL;
        }
        tabulatedFunctions[i] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
    }
}

ReferenceExtendedCustomCVForce::~ReferenceExtendedCustomCVForce() {
    for (int i = 0; i < tabulatedFunctions.size(); i++)
        if (tabulatedFunctions[i] != NULL)
            delete tabulatedFunctions[i];
}

void ReferenceExtendedCustomCVForce::calculateIxn(ContextImpl& innerContext, vector<Vec3>& atomCoordinates,
                                          const map<string, double>& globalParameters, vector<Vec3>& forces,
                                          double* totalEnergy, map<string, double>& energyParamDerivs) {
    // Compute the collective variables, and their derivatives with respect to particle positions.

    int numCVs = variableNames.size();
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(innerContext.getPlatformData());
    vector<Vec3>& innerForces = *((vector<Vec3>*) data->forces);
    map<string, double>& innerDerivs = *((map<string, double>*) data->energyParameterDerivatives);
    vector<vector<Vec3> > cvForces;
    vector<map<string, double> > cvDerivs;
    for (int i = 0; i < numCVs; i++) {
        cvValues[i] = innerContext.calcForcesAndEnergy(true, true, 1<<i);
        cvForces.push_back(innerForces);
        cvDerivs.push_back(innerDerivs);
    }

    // Compute the energy and forces.

    for (int i = 0; i < globalParameterNames.size(); i++)
        globalValues[i] = globalParameters.at(globalParameterNames[i]);
    int numParticles = atomCoordinates.size();
    if (totalEnergy != NULL)
        *totalEnergy += energyExpression.evaluate();
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate();
        for (int j = 0; j < numParticles; j++)
            forces[j] += cvForces[i][j]*dEdV;
    }

    // Compute the energy parameter derivatives.

    if (paramDerivExpressions.size() > 0) {
        for (int i = 0; i < paramDerivExpressions.size(); i++)
            energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate();
        for (int i = 0; i < numCVs; i++) {
            double dEdV = variableDerivExpressions[i].evaluate();
            for (auto& deriv : cvDerivs[i])
                energyParamDerivs[deriv.first] += dEdV*deriv.second;
        }
    }
}
