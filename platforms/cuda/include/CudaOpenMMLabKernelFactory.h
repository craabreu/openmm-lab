#ifndef OPENMM_CUDAOPENMM_LABKERNELFACTORY_H_
#define OPENMM_CUDAOPENMM_LABKERNELFACTORY_H_

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

using namespace OpenMM;

namespace OpenMMLab {

/**
 * This KernelFactory creates kernels for the CUDA implementation of the OpenMMLab plugin.
 */

class CudaOpenMMLabKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAOPENMM_LABKERNELFACTORY_H_*/
