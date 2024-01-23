/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2024 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */


#include "internal/ConcertedRMSDForceImpl.h"

#include "openmm/internal/CustomCPPForceImpl.h"
#include <vector>

using namespace OpenMMLab;
using namespace OpenMM;
using namespace std;

ConcertedRMSDForceImpl::ConcertedRMSDForceImpl(const ConcertedRMSDForce& owner) : CustomCPPForceImpl(owner), owner(owner) {
}

double ConcertedRMSDForceImpl::computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces) {
    return 0.0;
}
