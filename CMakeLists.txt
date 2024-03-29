#---------------------------------------------------
# OpenMM OpenMMLab Plugin
#----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

PROJECT(OpenMMOpenMMLab VERSION 0.0.0)

IF(DEFINED ENV{CONDA_PREFIX})
    SET(OPENMM_DIR $ENV{CONDA_PREFIX} CACHE PATH "Where OpenMM is installed")
ENDIF(DEFINED ENV{CONDA_PREFIX})

INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# Specify the C++ version we are building for.
SET (CMAKE_CXX_STANDARD 11)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

# Determine what architecture we are compiling for.
INCLUDE(TargetArch)
target_architecture(TARGET_ARCH)
if ("${TARGET_ARCH}" MATCHES "x86_64|i386")
    set(X86 ON)
endif()
if ("${TARGET_ARCH}" MATCHES "arm")
    set(ARM ON)
    if ("${TARGET_ARCH}" MATCHES "armv8")
        add_definitions(-D__ARM64__=1)
    else()
        add_definitions(-D__ARM__=1)
    endif()
endif()
if ("${TARGET_ARCH}" MATCHES "ppc")
    set(PPC ON)
    add_definitions(-D__PPC__=1)
endif()

# Set flags for linking on mac
IF(APPLE)
    SET(CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-msse2 -stdlib=libc++")
ENDIF(APPLE)

IF((X86 OR ARM) AND NOT (WIN32 AND OPENMM_BUILD_STATIC_LIB))
    SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DLEPTON_USE_JIT")
ENDIF()

# Select where to install
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    IF(WIN32)
        SET(CMAKE_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install the plugin" FORCE)
    ELSE(WIN32)
        SET(CMAKE_INSTALL_PREFIX $ENV{CONDA_PREFIX} CACHE PATH "Where to install the plugin" FORCE)
    ENDIF(WIN32)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(PLUGIN_SOURCE_SUBDIRS openmmapi serialization)

# Set the library name
SET(PLUGIN_LIBRARY_NAME OpenMMLab)
SET(SHARED_OPENMM_LAB_TARGET ${PLUGIN_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Create the library.

ADD_LIBRARY(${SHARED_OPENMM_LAB_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_OPENMM_LAB_TARGET}
    PROPERTIES COMPILE_FLAGS "-DPLUGIN_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_OPENMM_LAB_TARGET} OpenMM)
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_OPENMM_LAB_TARGET})

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

# Enable testing

ENABLE_TESTING()
ADD_SUBDIRECTORY(serialization/tests)
INCLUDE(CTest)

# Build the implementations for different platforms

ADD_SUBDIRECTORY(platforms/reference)
ADD_SUBDIRECTORY(platforms/common)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

FIND_PACKAGE(OPENCL QUIET)
IF(OPENCL_FOUND)
    MESSAGE("Found OpenCL: TRUE")
    SET(PLUGIN_BUILD_OPENCL_LIB ON CACHE BOOL "Build implementation for OpenCL")
ELSE(OPENCL_FOUND)
    MESSAGE("Found OpenCL: FALSE")
    SET(PLUGIN_BUILD_OPENCL_LIB OFF CACHE BOOL "Build implementation for OpenCL")
ENDIF(OPENCL_FOUND)
IF(PLUGIN_BUILD_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(PLUGIN_BUILD_OPENCL_LIB)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    MESSAGE("Found CUDA: TRUE")
    SET(PLUGIN_BUILD_CUDA_LIB ON CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    MESSAGE("Found CUDA: FALSE")
    SET(PLUGIN_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(PLUGIN_BUILD_CUDA_LIB)
    GET_FILENAME_COMPONENT(CUDA_LIB_DIR "${CUDA_cufft_LIBRARY}" PATH)
    SET(CUDA_nvrtc_LIBRARY "${CUDA_LIB_DIR}/libnvrtc.so")
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(PLUGIN_BUILD_CUDA_LIB)

# Fetch VkFFT

IF(OPENCL_FOUND OR CUDA_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        VkFFT
        GIT_REPOSITORY https://github.com/DTolm/VkFFT
        GIT_TAG        v1.3.2
        CMAKE_ARGS     -DVKFFT_BACKEND=3
    )
    FetchContent_Populate(vkFFT)
ENDIF(OPENCL_FOUND OR CUDA_FOUND)

# Build the Python API

FIND_PROGRAM(DOXYGEN_EXECUTABLE doxygen)
FIND_PROGRAM(PYTHON_EXECUTABLE python REQUIRED)
FIND_PROGRAM(SWIG_EXECUTABLE swig REQUIRED)
SET(PLUGIN_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
IF(PLUGIN_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(PLUGIN_BUILD_PYTHON_WRAPPERS)

# Build the documentation

ADD_SUBDIRECTORY(docs)

# Find required packages

SET(CMAKE_THREAD_PREFER_PTHREAD ON CACHE BOOL "Prefer pthread lib")
SET(THREADS_PREFER_PTHREAD_FLAG ON CACHE BOOL "Prefer pthread flag")
FIND_PACKAGE(Threads REQUIRED)
