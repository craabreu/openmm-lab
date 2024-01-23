#ifndef OPENMM_CONCERTEDRMSDFORCEIMPL_H_
#define OPENMM_CONCERTEDRMSDFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2024 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "ConcertedRMSDForce.h"

#include "openmm/internal/CustomCPPForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include <vector>

using namespace OpenMM;
using namespace std;

namespace OpenMMLab {

/**
 * This is the internal implementation of ConcertedRMSDForce.
 */

class ConcertedRMSDForceImpl : public CustomCPPForceImpl {
public:
    ConcertedRMSDForceImpl(const ConcertedRMSDForce& owner) : CustomCPPForceImpl(owner), owner(owner) {}
    void initialize(ContextImpl& context);
    double computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces);
    const ConcertedRMSDForce& getOwner() const {
        return owner;
    }
private:
    const ConcertedRMSDForce& owner;
    int numParticles;
    vector<int> particles;
    vector<Vec3> referencePos;
};

} // namespace OpenMMLab

#endif /*OPENMM_CONCERTEDRMSDFORCEIMPL_H_*/
