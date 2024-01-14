#ifndef COMMON_OPENMM_LAB_KERNELS_H_
#define COMMON_OPENMM_LAB_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "OpenMMLabKernels.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeKernel.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionProgram.h"

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This kernel is invoked by ExtendedCustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcExtendedCustomCVForceKernel : public CalcExtendedCustomCVForceKernel {
public:
    CommonCalcExtendedCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CalcExtendedCustomCVForceKernel(name, platform),
            cc(cc), hasInitializedListeners(false) {
    }
    ~CommonCalcExtendedCustomCVForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the ExtendedCustomCVForce this kernel will be used for
     * @param innerContext   the context created by the ExtendedCustomCVForce for computing collective variables
     */
    void initialize(const System& system, const ExtendedCustomCVForce& force, ContextImpl& innerContext);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the ExtendedCustomCVForce for computing collective variables
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy);
    /**
     * Copy state information to the inner context.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the ExtendedCustomCVForce for computing collective variables
     */
    void copyState(ContextImpl& context, ContextImpl& innerContext);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ExtendedCustomCVForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const ExtendedCustomCVForce& force);
    /**
     * Get the ComputeContext corresponding to the inner Context.
     */
    virtual ComputeContext& getInnerComputeContext(ContextImpl& innerContext) = 0;
private:
    class ForceInfo;
    class ReorderListener;
    class TabulatedFunctionWrapper;
    ComputeContext& cc;
    bool hasInitializedListeners;
    Lepton::CompiledExpression energyExpression;
    std::vector<std::string> variableNames, paramDerivNames, globalParameterNames;
    std::vector<Lepton::CompiledExpression> variableDerivExpressions;
    std::vector<Lepton::CompiledExpression> paramDerivExpressions;
    std::vector<ComputeArray> cvForces;
    std::vector<double> globalValues, cvValues;
    std::vector<Lepton::CustomFunction*> tabulatedFunctions;
    ComputeArray invAtomOrder;
    ComputeArray innerInvAtomOrder;
    ComputeKernel copyStateKernel, copyForcesKernel, addForcesKernel;
};

} // namespace OpenMMLab

#endif /*COMMON_OPENMM_LAB_KERNELS_H_*/
