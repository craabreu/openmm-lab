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
#include <vector>

using namespace OpenMM;
using namespace std;

namespace OpenMMLab {

/**
 * This is the internal implementation of ConcertedRMSDForce.
 */

class ConcertedRMSDForceImpl : public CustomCPPForceImpl {
public:
    ConcertedRMSDForceImpl(const ConcertedRMSDForce& owner);
    double computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces);
    const ConcertedRMSDForce& getOwner() const {
        return owner;
    }
private:
    const ConcertedRMSDForce& owner;
};

} // namespace OpenMMLab

#endif /*OPENMM_CONCERTEDRMSDFORCEIMPL_H_*/
