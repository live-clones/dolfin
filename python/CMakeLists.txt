cmake_minimum_required(VERSION 3.5.0)

PROJECT(dolfin_testYYY)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})
add_subdirectory(pybind11)

# Create the binding library
pybind11_add_module(cpp SHARED
  src/dolfin.cpp
  src/experimental/experimental.cpp
  src/function/function.cpp
  src/generation/generation.cpp
  src/geometry/geometry.cpp
  src/mesh/mesh.cpp
  src/io/io.cpp
  )

# Add DOLFIN libraries and other config
target_link_libraries(cpp PRIVATE pybind11::module dolfin)

#set_target_properties(dolfin_test PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
#  SUFFIX "${PYTHON_MODULE_EXTENSION}")
