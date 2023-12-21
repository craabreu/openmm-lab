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
#include "openmm/reference/ReferencePlatform.h"
#include <map>
#include <string>

extern "C" OPENMM_EXPORT void registerOpenMMLabReferenceKernelFactories();

OpenMM::ReferencePlatform platform;
std::map<std::string, std::string> properties;

void initializeTests(int argc, char* argv[]) {
    registerOpenMMLabReferenceKernelFactories();
    platform = dynamic_cast<OpenMM::ReferencePlatform&>(OpenMM::Platform::getPlatformByName("Reference"));
}