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

ReferenceExtendedCustomCVForce::ReferenceExtendedCustomCVForce(const ExtendedCustomCVForce& force) {
    for (int i = 0; i < force.getNumCollectiveVariables(); i++)
        variableNames.push_back(force.getCollectiveVariableName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        paramDerivNames.push_back(force.getEnergyParameterDerivativeName(i));

    // Create custom functions for the tabulated functions.

    map<string, CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Create the expressions.

    ParsedExpression energyExpr = Parser::parse(force.getEnergyFunction(), functions);
    energyExpression = energyExpr.createProgram();
    variableDerivExpressions.clear();
    for (auto& name : variableNames)
        variableDerivExpressions.push_back(energyExpr.differentiate(name).optimize().createProgram());
    paramDerivExpressions.clear();
    for (auto& name : paramDerivNames)
        paramDerivExpressions.push_back(energyExpr.differentiate(name).optimize().createProgram());

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

static void replaceFunctionsInExpression(map<string, CustomFunction*>& functions, ExpressionProgram& expression) {
    for (int i = 0; i < expression.getNumOperations(); i++) {
        if (expression.getOperation(i).getId() == Operation::CUSTOM) {
            const Operation::Custom& op = dynamic_cast<const Operation::Custom&>(expression.getOperation(i));
            expression.setOperation(i, new Operation::Custom(op.getName(), functions[op.getName()]->clone(), op.getDerivOrder()));
        }
    }
}

void ReferenceExtendedCustomCVForce::updateTabulatedFunctions(const ExtendedCustomCVForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Replace tabulated functions in the expressions.

    replaceFunctionsInExpression(functions, energyExpression);
    for (auto& expression : variableDerivExpressions)
        replaceFunctionsInExpression(functions, expression);
    for (auto& expression : paramDerivExpressions)
        replaceFunctionsInExpression(functions, expression);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

ReferenceExtendedCustomCVForce::~ReferenceExtendedCustomCVForce() {
}

void ReferenceExtendedCustomCVForce::calculateIxn(ContextImpl& innerContext, vector<Vec3>& atomCoordinates,
                                          const map<string, double>& globalParameters, vector<Vec3>& forces,
                                          double* totalEnergy, map<string, double>& energyParamDerivs) const {
    // Compute the collective variables, and their derivatives with respect to particle positions.

    int numCVs = variableNames.size();
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(innerContext.getPlatformData());
    vector<Vec3>& innerForces = *((vector<Vec3>*) data->forces);
    map<string, double>& innerDerivs = *((map<string, double>*) data->energyParameterDerivatives);
    vector<double> cvValues;
    vector<vector<Vec3> > cvForces;
    vector<map<string, double> > cvDerivs;
    for (int i = 0; i < numCVs; i++) {
        cvValues.push_back(innerContext.calcForcesAndEnergy(true, true, 1<<i));
        cvForces.push_back(innerForces);
        cvDerivs.push_back(innerDerivs);
    }

    // Compute the energy and forces.

    int numParticles = atomCoordinates.size();
    map<string, double> variables = globalParameters;
    for (int i = 0; i < numCVs; i++)
        variables[variableNames[i]] = cvValues[i];
    if (totalEnergy != NULL)
        *totalEnergy += energyExpression.evaluate(variables);
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate(variables);
        for (int j = 0; j < numParticles; j++)
            forces[j] += cvForces[i][j]*dEdV;
    }

    // Compute the energy parameter derivatives.

    for (int i = 0; i < paramDerivExpressions.size(); i++)
        energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate(variables);
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate(variables);
        for (auto& deriv : cvDerivs[i])
            energyParamDerivs[deriv.first] += dEdV*deriv.second;
    }
}
