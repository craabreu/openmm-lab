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
#include <windows.h>
#include <sstream>
#else
#include <dlfcn.h>
#include <dirent.h>
#include <cstdlib>
#endif

#include "SlicedNonbondedForce.h"
#include "SlicedNonbondedForceProxy.h"
#include "ConcertedRMSDForce.h"
#include "ConcertedRMSDForceProxy.h"
#include "openmm/serialization/SerializationProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" OPENMM_EXPORT_OPENMM_LAB void registerOpenMMLabSerializationProxies();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerOpenMMLabSerializationProxies();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) registerOpenMMLabSerializationProxies();
#endif

using namespace OpenMMLab;
using namespace OpenMM;

extern "C" OPENMM_EXPORT_OPENMM_LAB void registerOpenMMLabSerializationProxies() {
    SerializationProxy::registerProxy(typeid(SlicedNonbondedForce), new SlicedNonbondedForceProxy());
    SerializationProxy::registerProxy(typeid(ConcertedRMSDForce), new ConcertedRMSDForceProxy());
}
