#ifndef OPENMM_OPENCLOPENMM_LABKERNELFACTORY_H_
#define OPENMM_OPENCLOPENMM_LABKERNELFACTORY_H_

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

namespace OpenMMLab {

/**
 * This KernelFactory creates kernels for the OpenCL implementation of the OpenMMLab plugin.
 */

class OpenCLOpenMMLabKernelFactory : public OpenMM::KernelFactory {
public:
    OpenMM::KernelImpl* createKernelImpl(std::string name, const OpenMM::Platform& platform, OpenMM::ContextImpl& context) const;
};

} // namespace OpenMMLab

#endif /*OPENMM_OPENCLOPENMM_LABKERNELFACTORY_H_*/
