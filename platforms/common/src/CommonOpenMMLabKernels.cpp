/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "CommonOpenMMLabKernels.h"
#include "CommonOpenMMLabKernelSources.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ContextImpl.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "openmm/reference/ReferenceTabulatedFunction.h"

using namespace OpenMMLab;
using namespace OpenMM;
using namespace Lepton;
using namespace std;


class CommonCalcExtendedCustomCVForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(ComputeForceInfo& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        return force.areParticlesIdentical(particle1, particle2);
    }
    int getNumParticleGroups() {
        return force.getNumParticleGroups();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        force.getParticlesInGroup(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return force.areGroupsIdentical(group1, group2);
    }
private:
    ComputeForceInfo& force;
};

class CommonCalcExtendedCustomCVForceKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(ComputeContext& cc, ArrayInterface& invAtomOrder) : cc(cc), invAtomOrder(invAtomOrder) {
    }
    void execute() {
        vector<int> invOrder(cc.getPaddedNumAtoms());
        const vector<int>& order = cc.getAtomIndex();
        for (int i = 0; i < order.size(); i++)
            invOrder[order[i]] = i;
        invAtomOrder.upload(invOrder);
    }
private:
    ComputeContext& cc;
    ArrayInterface& invAtomOrder;
};

// This class allows us to update tabulated functions without having to recompile expressions
// that use them.
class CommonCalcExtendedCustomCVForceKernel::TabulatedFunctionWrapper : public CustomFunction {
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

void CommonCalcExtendedCustomCVForceKernel::initialize(const System& system, const ExtendedCustomCVForce& force, ContextImpl& innerContext) {
    ContextSelector selector(cc);
    int numCVs = force.getNumCollectiveVariables();
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < numCVs; i++)
        variableNames.push_back(force.getCollectiveVariableName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string name = force.getEnergyParameterDerivativeName(i);
        paramDerivNames.push_back(name);
        cc.addEnergyParameterDerivative(name);
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions(), NULL);
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        tabulatedFunctions[i] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
        functions[force.getTabulatedFunctionName(i)] = new TabulatedFunctionWrapper(tabulatedFunctions, i);
    }

    // Create the expressions.

    Lepton::ParsedExpression energyExpr = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
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

    // Copy parameter derivatives from the inner context.

    ComputeContext& cc2 = getInnerComputeContext(innerContext);
    for (auto& param : cc2.getEnergyParamDerivNames())
        cc.addEnergyParameterDerivative(param);

    // Create arrays for storing information.

    cvForces.resize(numCVs);
    for (int i = 0; i < numCVs; i++)
        cvForces[i].initialize<long long>(cc, 3*cc.getPaddedNumAtoms(), "cvForce");
    invAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "invAtomOrder");
    innerInvAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "innerInvAtomOrder");

    // Create the kernels.

    stringstream args, add;
    for (int i = 0; i < numCVs; i++) {
        args << ", GLOBAL mm_long * RESTRICT force" << i << ", real dEdV" << i;
        add << "forces[i] += (mm_long) (force" << i << "[i]*dEdV" << i << ");\n";
    }
    map<string, string> replacements;
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    replacements["ADD_FORCES"] = add.str();
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonOpenMMLabKernelSources::customCVForce, replacements));
    copyStateKernel = program->createKernel("copyState");
    copyStateKernel->addArg(cc.getPosq());
    copyStateKernel->addArg(cc2.getPosq());
    if (cc.getUseMixedPrecision()) {
        copyStateKernel->addArg(cc.getPosqCorrection());
        copyStateKernel->addArg(cc2.getPosqCorrection());
    }
    copyStateKernel->addArg(cc.getVelm());
    copyStateKernel->addArg(cc2.getVelm());
    copyStateKernel->addArg(cc.getAtomIndexArray());
    copyStateKernel->addArg(innerInvAtomOrder);
    copyStateKernel->addArg(cc.getNumAtoms());
    copyForcesKernel = program->createKernel("copyForces");
    copyForcesKernel->addArg();
    copyForcesKernel->addArg(invAtomOrder);
    copyForcesKernel->addArg(cc2.getLongForceBuffer());
    copyForcesKernel->addArg(cc2.getAtomIndexArray());
    copyForcesKernel->addArg(cc.getNumAtoms());
    copyForcesKernel->addArg(cc.getPaddedNumAtoms());
    addForcesKernel = program->createKernel("addForces");
    addForcesKernel->addArg(cc.getLongForceBuffer());
    addForcesKernel->addArg((int) cc.getLongForceBuffer().getSize());
    for (int i = 0; i < numCVs; i++) {
        addForcesKernel->addArg();
        addForcesKernel->addArg();
    }

    // This context needs to respect all forces in the inner context when reordering atoms.

    for (auto* info : cc2.getForceInfos())
        cc.addForce(new ForceInfo(*info));
}

CommonCalcExtendedCustomCVForceKernel::~CommonCalcExtendedCustomCVForceKernel() {
    for (int i = 0; i < tabulatedFunctions.size(); i++)
        if (tabulatedFunctions[i] != NULL)
            delete tabulatedFunctions[i];
}

double CommonCalcExtendedCustomCVForceKernel::execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy) {
    copyState(context, innerContext);
    int numCVs = variableNames.size();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    vector<map<string, double> > cvDerivs(numCVs);
    for (int i = 0; i < numCVs; i++) {
        cvValues[i] = innerContext.calcForcesAndEnergy(true, true, 1<<i);
        ContextSelector selector(cc);
        copyForcesKernel->setArg(0, cvForces[i]);
        copyForcesKernel->execute(numAtoms);
        innerContext.getEnergyParameterDerivatives(cvDerivs[i]);
    }

    // Compute the energy and forces.

    ContextSelector selector(cc);
    for (int i = 0; i < globalParameterNames.size(); i++)
        globalValues[i] = context.getParameter(globalParameterNames[i]);
    double energy = energyExpression.evaluate();
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate();
        addForcesKernel->setArg(2*i+2, cvForces[i]);
        if (cc.getUseDoublePrecision())
            addForcesKernel->setArg(2*i+3, dEdV);
        else
            addForcesKernel->setArg(2*i+3, (float) dEdV);
    }
    addForcesKernel->execute(numAtoms);

    // Compute the energy parameter derivatives.

    if (paramDerivExpressions.size() > 0) {
        map<string, double>& energyParamDerivs = cc.getEnergyParamDerivWorkspace();
        for (int i = 0; i < paramDerivExpressions.size(); i++)
            energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate();
        for (int i = 0; i < numCVs; i++) {
            double dEdV = variableDerivExpressions[i].evaluate();
            for (auto& deriv : cvDerivs[i])
                energyParamDerivs[deriv.first] += dEdV*deriv.second;
        }
    }
    return energy;
}

void CommonCalcExtendedCustomCVForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext) {
    ContextSelector selector(cc);
    int numAtoms = cc.getNumAtoms();
    ComputeContext& cc2 = getInnerComputeContext(innerContext);
    if (!hasInitializedListeners) {
        hasInitializedListeners = true;

        // Initialize the listeners.

        ReorderListener* listener1 = new ReorderListener(cc, invAtomOrder);
        ReorderListener* listener2 = new ReorderListener(cc2, innerInvAtomOrder);
        cc.addReorderListener(listener1);
        cc2.addReorderListener(listener2);
        listener1->execute();
        listener2->execute();
    }
    copyStateKernel->execute(numAtoms);
    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext.setPeriodicBoxVectors(a, b, c);
    innerContext.setTime(context.getTime());
    map<string, double> innerParameters = innerContext.getParameters();
    for (auto& param : innerParameters)
        innerContext.setParameter(param.first, context.getParameter(param.first));
}

void CommonCalcExtendedCustomCVForceKernel::copyParametersToContext(ContextImpl& context, const ExtendedCustomCVForce& force) {
    // Create custom functions for the tabulated functions.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        if (tabulatedFunctions[i] != NULL) {
            delete tabulatedFunctions[i];
            tabulatedFunctions[i] = NULL;
        }
        tabulatedFunctions[i] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
    }
}
