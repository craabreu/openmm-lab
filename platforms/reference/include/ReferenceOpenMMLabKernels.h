#ifndef REFERENCE_OPENMM_LAB_KERNELS_H_
#define REFERENCE_OPENMM_LAB_KERNELS_H_

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
#include "ReferenceExtendedCustomCVForce.h"
#include "openmm/Platform.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include <vector>
#include <array>
#include <map>

using namespace std;

namespace OpenMMLab {

/**
 * This kernel is invoked by SlicedNonbondedForce to calculate the forces acting on the system.
 */
class ReferenceCalcSlicedNonbondedForceKernel : public CalcSlicedNonbondedForceKernel {
public:
    ReferenceCalcSlicedNonbondedForceKernel(string name, const Platform& platform) : CalcSlicedNonbondedForceKernel(name, platform) {
    }
    ~ReferenceCalcSlicedNonbondedForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the SlicedNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const SlicedNonbondedForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param includeReciprocal  true if reciprocal space interactions should be included
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the SlicedNonbondedForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const SlicedNonbondedForce& force);
    /**
     * Get the parameters being used for PME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Get the dispersion parameters being used for the dispersion term in LJPME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
    static const int Coul = 0;
    static const int vdW = 1;
    class ScalingParameterInfo;
    void computeParameters(ContextImpl& context);
    int numParticles, num14;
    vector<vector<int>>bonded14IndexArray;
    vector<vector<double>> particleParamArray, bonded14ParamArray;
    vector<int> bonded14SliceArray;
    vector<array<double, 3>> baseParticleParams, baseExceptionParams;
    map<pair<string, int>, array<double, 3>> particleParamOffsets, exceptionParamOffsets;
    double nonbondedCutoff, switchingDistance, rfDielectric, ewaldAlpha, ewaldDispersionAlpha;
    vector<double> dispersionCoefficients;
    int kmax[3], gridSize[3], dispersionGridSize[3];
    bool useSwitchingFunction, exceptionsArePeriodic;
    vector<set<int>> exclusions;
    NonbondedMethod nonbondedMethod;
    NeighborList* neighborList;

    int numSubsets, numSlices;
    vector<int> subsets;
    vector<vector<double>> sliceLambdas;
    vector<vector<ScalingParameterInfo>> sliceScalingParams;
};

class ReferenceCalcSlicedNonbondedForceKernel::ScalingParameterInfo {
public:
    string name;
    bool hasDerivative;
    ScalingParameterInfo() : name(""), hasDerivative(false) {}
    ScalingParameterInfo(string name, bool hasDerivative) : name(name), hasDerivative(hasDerivative) {}
};

/**
 * This kernel is invoked by ExtendedCustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcExtendedCustomCVForceKernel : public CalcExtendedCustomCVForceKernel {
public:
    ReferenceCalcExtendedCustomCVForceKernel(std::string name, const Platform& platform) : CalcExtendedCustomCVForceKernel(name, platform), ixn(NULL) {
    }
    ~ReferenceCalcExtendedCustomCVForceKernel();
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
private:
    ReferenceExtendedCustomCVForce* ixn;
    std::vector<std::string> globalParameterNames, energyParamDerivNames;
};

} // namespace OpenMMLab

#endif /*REFERENCE_OPENMM_LAB_KERNELS_H_*/
