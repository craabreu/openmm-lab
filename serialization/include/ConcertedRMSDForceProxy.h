#ifndef OPENMM_CONCERTEDRMSDFORCE_PROXY_H_
#define OPENMM_CONCERTEDRMSDFORCE_PROXY_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2024 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "internal/windowsExportOpenMMLab.h"

#include "openmm/serialization/SerializationProxy.h"

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This is a proxy for serializing ConcertedRMSDForce objects.
 */

class OPENMM_EXPORT_OPENMM_LAB ConcertedRMSDForceProxy : public SerializationProxy {
public:
    ConcertedRMSDForceProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

} // namespace OpenMMLab

#endif /*OPENMM_CONCERTEDRMSDFORCE_PROXY_H_*/
