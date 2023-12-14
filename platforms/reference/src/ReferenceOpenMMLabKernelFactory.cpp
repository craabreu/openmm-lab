/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceOpenMMLabKernelFactory.h"
#include "ReferenceOpenMMLabKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMMLab;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
            ReferenceOpenMMLabKernelFactory* factory = new ReferenceOpenMMLabKernelFactory();
            platform.registerKernelFactory(CalcSlicedNonbondedForceKernel::Name(), factory);
            platform.registerKernelFactory(CalcExtendedCustomCVForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerOpenMMLabReferenceKernelFactories() {
    registerKernelFactories();
}

KernelImpl* ReferenceOpenMMLabKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    if (name == CalcSlicedNonbondedForceKernel::Name())
        return new ReferenceCalcSlicedNonbondedForceKernel(name, platform);
    if (name == CalcExtendedCustomCVForceKernel::Name())
        return new ReferenceCalcExtendedCustomCVForceKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
