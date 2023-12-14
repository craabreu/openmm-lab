/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include <exception>

#include "OpenCLOpenMMLabKernelFactory.h"
#include "OpenCLOpenMMLabKernels.h"
#include "OpenCLParallelOpenMMLabKernels.h"
#include "CommonOpenMMLabKernels.h"
#include "openmm/opencl/OpenCLContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMMLab;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("OpenCL");
        OpenCLOpenMMLabKernelFactory* factory = new OpenCLOpenMMLabKernelFactory();
        platform.registerKernelFactory(CalcSlicedNonbondedForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcExtendedCustomCVForceKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerOpenMMLabOpenCLKernelFactories() {
    try {
        Platform::getPlatformByName("OpenCL");
    }
    catch (...) {
        Platform::registerPlatform(new OpenCLPlatform());
    }
    registerKernelFactories();
}

KernelImpl* OpenCLOpenMMLabKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    OpenCLPlatform::PlatformData& data = *static_cast<OpenCLPlatform::PlatformData*>(context.getPlatformData());
    if (data.contexts.size() > 1) {
        if (name == CalcSlicedNonbondedForceKernel::Name())
            return new OpenCLParallelCalcSlicedNonbondedForceKernel(name, platform, data, context.getSystem());
        throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
    }
    OpenCLContext& cl = *data.contexts[0];
    if (name == CalcSlicedNonbondedForceKernel::Name())
        return new OpenCLCalcSlicedNonbondedForceKernel(name, platform, cl, context.getSystem());
    if (name == CalcExtendedCustomCVForceKernel::Name())
        return new OpenCLCalcExtendedCustomCVForceKernel(name, platform, cl);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
