/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/cuda/CudaPlatform.h"

extern "C" OPENMM_EXPORT void registerOpenMMLabCudaKernelFactories();

OpenMM::CudaPlatform platform;

void initializeTests(int argc, char* argv[]) {
    registerOpenMMLabCudaKernelFactories();
    platform = dynamic_cast<OpenMM::CudaPlatform&>(OpenMM::Platform::getPlatformByName("CUDA"));
    if (argc > 1)
        platform.setPropertyDefaultValue("Precision", std::string(argv[1]));
    if (argc > 2)
        platform.setPropertyDefaultValue("DeviceIndex", std::string(argv[2]));
}