OpenMM Laboratory Plugin
========================

[![Linux](https://github.com/craabreu/openmm-lab/actions/workflows/Linux.yml/badge.svg)](https://github.com/craabreu/openmm-lab/actions/workflows/Linux.yml)
[![MacOS](https://github.com/craabreu/openmm-lab/actions/workflows/MacOS.yml/badge.svg)](https://github.com/craabreu/openmm-lab/actions/workflows/MacOS.yml)
[![Doc](https://github.com/craabreu/openmm-lab/actions/workflows/Doc.yml/badge.svg)](https://github.com/craabreu/openmm-lab/actions/workflows/Doc.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This [OpenMM] plugin is a laboratory for low-level code implementation for OpenMM.

Documentation
=============

Documentation for this plugin is available at [Github Pages](https://craabreu.github.io/openmm-lab/).
It includes the Python API and the theory for slicing lattice-sum energy contributions.

Installing from Source
======================

This project uses [CMake] for its build system.  To build it, follow these steps:

1. Create a directory in which to build the plugin.

2. Run the CMake GUI or ccmake, specifying your new directory as the build directory and the top
level directory of this project as the source directory.

3. Press "Configure".

4. Set OPENMM_DIR to point to the directory where OpenMM is installed.  This is needed to locate
the OpenMM header files and libraries.

5. Set CMAKE_INSTALL_PREFIX to the directory where the plugin should be installed.  Usually,
this will be the same as OPENMM_DIR, so the plugin will be added to your OpenMM installation.

6. If you plan to build the OpenCL platform, make sure that OPENCL_INCLUDE_DIR and
OPENCL_LIBRARY are set correctly, and that PLUGIN_BUILD_OPENCL_LIB is selected.

7. If you plan to build the CUDA platform, make sure that CUDA_TOOLKIT_ROOT_DIR is set correctly
and that PLUGIN_BUILD_CUDA_LIB is selected.

8. Press "Configure" again if necessary, then press "Generate".

9. Use the build system you selected to build and install the plugin.  For example, if you
selected Unix Makefiles, type `make install`.

Python Wrapper and API
======================

As [OpenMM], this project uses [SWIG] to generate its Python API.  SWIG takes an "interface
file", which is essentially a C++ header file with some extra annotations added, as its input.
It then generates a Python extension module exposing the C++ API in Python.

To build and install the Python API, build the `PythonInstall` target, for example by typing
`make PythonInstall` (if you are installing into the system Python, you may need to use sudo).
Once you do that, you can use the plugin from your Python scripts:

```py
    import openmm as mm
    import openmmlab as nbs
    system = mm.System()
    force = nbs.SlicedNonbondedForce(2)
    system.addForce(force)
```

Test Cases
==========

To run the C++ test cases, build the "test" target, for example by typing `make test`.

To run the Python test cases, build the "PythonTest" target by typing `make PythonTest`.


[CMake]:                http://www.cmake.org
[NonbondedForce]:       http://docs.openmm.org/latest/api-python/generated/openmm.openmm.NonbondedForce.html
[Context]:              http://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html
[getState]:             http://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context.getState
[OpenMM]:               https://openmm.org
[SWIG]:                 http://www.swig.org