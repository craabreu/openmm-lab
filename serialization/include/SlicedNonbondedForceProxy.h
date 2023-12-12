#ifndef OPENMM_SLICEDNONBONDEDFORCE_PROXY_H_
#define OPENMM_SLICEDNONBONDEDFORCE_PROXY_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "internal/windowsExportOpenMMLab.h"
#include "openmm/serialization/SerializationProxy.h"

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This is a proxy for serializing SlicedNonbondedForce objects.
 */

class OPENMM_EXPORT_OPENMM_LAB SlicedNonbondedForceProxy : public SerializationProxy {
public:
    SlicedNonbondedForceProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

} // namespace OpenMM

#endif /*OPENMM_SLICEDNONBONDEDFORCE_PROXY_H_*/
