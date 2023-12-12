#ifndef OPENMM_SLICEDNONBONDEDFORCEIMPL_H_
#define OPENMM_SLICEDNONBONDEDFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "SlicedNonbondedForce.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This is the internal implementation of SlicedNonbondedForce.
 */

class OPENMM_EXPORT_OPENMM_LAB SlicedNonbondedForceImpl : public NonbondedForceImpl {
public:
    SlicedNonbondedForceImpl(const SlicedNonbondedForce& owner);
    ~SlicedNonbondedForceImpl();
    void initialize(ContextImpl& context);
    const SlicedNonbondedForce& getOwner() const {
        return owner;
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    static double calcDispersionCorrection(const System& system, const SlicedNonbondedForce& force);
    static vector<double> calcDispersionCorrections(const System& system, const SlicedNonbondedForce& force);
private:
    static double evalIntegral(double r, double rs, double rc, double sigma);
    const SlicedNonbondedForce& owner;
    Kernel kernel;
};

} // namespace OpenMMLab

#endif /*OPENMM_SLICEDNONBONDEDFORCEIMPL_H_*/
