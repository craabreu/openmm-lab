from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
openmmlab_header_dir = '@PLUGIN_HEADER_DIR@'
openmmlab_library_dir = '@PLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

os.environ['CC'] = '@CMAKE_C_COMPILER@'
os.environ['CXX'] = '@CMAKE_CXX_COMPILER@'

extension = Extension(
    name='_openmmlab',
    sources=['OpenMMLabWrapper.cpp'],
    libraries=['OpenMM', 'OpenMMLab'],
    include_dirs=[os.path.join(openmm_dir, 'include'), openmmlab_header_dir],
    library_dirs=[os.path.join(openmm_dir, 'lib'), openmmlab_library_dir],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
)

setup(
    name='openmmlab',
    version='@CMAKE_PROJECT_VERSION@',
    py_modules=['openmmlab'],
    ext_modules=[extension],
)
