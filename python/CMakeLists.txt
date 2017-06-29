cmake_minimum_required(VERSION 3.5.0)

PROJECT(dolfin_test)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")
#include(cmake-modules/pybind11Tools.cmake)
add_subdirectory(pybind11)
#find_package(pybind11 REQUIRED)


# Create the binding library
pybind11_add_module(dolfin_test SHARED
  src/dolfin.cpp
  src/generation/generation.cpp
  src/geometry/geometry.cpp
  src/mesh/mesh.cpp
  src/io/io.cpp
  )

target_link_libraries(dolfin_test PRIVATE pybind11::module dolfin)
set_target_properties(dolfin_test PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}")
