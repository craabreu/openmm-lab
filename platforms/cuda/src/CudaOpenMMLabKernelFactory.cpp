/* -------------------------------------------------------------------------- *
 *                              OpenMMOpenMMLab                                   *
 * -------------------------------------------------------------------------- */

#include <exception>

#include "CudaOpenMMLabKernelFactory.h"
#include "CudaOpenMMLabKernels.h"
#include "CudaParallelOpenMMLabKernels.h"
#include "CommonOpenMMLabKernels.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMMLab;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        CudaOpenMMLabKernelFactory* factory = new CudaOpenMMLabKernelFactory();
        platform.registerKernelFactory(CalcSlicedNonbondedForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcExtendedCustomCVForceKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerOpenMMLabCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}

KernelImpl* CudaOpenMMLabKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaPlatform::PlatformData& data = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    if (data.contexts.size() > 1) {
        // We are running in parallel on multiple devices, so we may want to create a parallel kernel.
        if (name == CalcSlicedNonbondedForceKernel::Name())
            return new CudaParallelCalcSlicedNonbondedForceKernel(name, platform, data, context.getSystem());
    }
    CudaContext& cu = *data.contexts[0];
    if (name == CalcSlicedNonbondedForceKernel::Name())
        return new CudaCalcSlicedNonbondedForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcExtendedCustomCVForceKernel::Name())
        return new CudaCalcExtendedCustomCVForceKernel(name, platform, cu);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
