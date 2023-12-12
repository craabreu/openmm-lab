#ifndef OPENMM_REFERENCEOPENMM_LABKERNELFACTORY_H_
#define OPENMM_REFERENCEOPENMM_LABKERNELFACTORY_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/KernelFactory.h"
#include <string.h>

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This KernelFactory creates kernels for the reference implementation of the OpenMMLab plugin.
 */

class ReferenceOpenMMLabKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMMLab

#endif /*OPENMM_REFERENCEOPENMM_LABKERNELFACTORY_H_*/
