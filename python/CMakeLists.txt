cmake_minimum_required(VERSION 3.5.0)

set(NAME example)
PROJECT(${NAME})

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

#find_package(PythonInterp 3)

# Get Python include path from Python interpretter
#execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
#  "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
#  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
#  RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)

# Get Python library path from interpreter
#execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
#  "import os, sys, inspect; sys.stdout.write(os.path.split(os.path.split(inspect.getfile(inspect))[0])[0])"
#  OUTPUT_VARIABLE _PYTHON_LIB_PATH
#  RESULT_VARIABLE _PYTHON_LIB_RESULT)

# Set include path, if returned by interpreter
#if ("${_PYTHON_INCLUDE_RESULT}" STREQUAL "0")
#  set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH})
#endif()

#find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
#INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

find_package(pybind11 REQUIRED)

#pybind11_add_module(dolfin_test.so
#  dolfin/dolfin.cpp
#  dolfin/generation/generation.cpp
#  dolfin/geometry/geometry.cpp
#  dolfin/mesh/mesh.cpp
#  dolfin/io/io.cpp
#  )

# Create the binding library
pybind11_add_module(dolfin_test SHARED
#add_library(dolfin_test MODULE
  dolfin/dolfin.cpp
  dolfin/generation/generation.cpp
  dolfin/geometry/geometry.cpp
  dolfin/mesh/mesh.cpp
  dolfin/io/io.cpp
  )

target_link_libraries(dolfin_test PRIVATE pybind11::module dolfin)
set_target_properties(dolfin_test PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}")
