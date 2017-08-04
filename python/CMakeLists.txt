cmake_minimum_required(VERSION 3.5.0)

PROJECT(dolfin_testYYY)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})
add_subdirectory(pybind11)


# Check for mpi4py
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#find_package(MPI4PY)

# Create the binding library
pybind11_add_module(cpp SHARED
  src/dolfin.cpp
  src/adaptivity/adaptivity.cpp
  src/ale/ale.cpp
  src/common/common.cpp
  src/fem/fem.cpp
  src/function/function.cpp
  src/generation/generation.cpp
  src/geometry/geometry.cpp
  src/graph/graph.cpp
  src/log/log.cpp
  src/math/math.cpp
  src/mesh/mesh.cpp
  src/multistage/multistage.cpp
  src/io/io.cpp
  src/la/la.cpp
  src/nls/nls.cpp
  src/parameter/parameter.cpp
  src/refinement/refinement.cpp
  )

if (MPI4PY_FOUND)
  target_include_directories(cpp PRIVATE ${MPI4PY_INCLUDE_DIR})
  add_definitions(-DHAS_MPI4PY)
endif()

# Add DOLFIN libraries and other config
target_link_libraries(cpp PRIVATE pybind11::module dolfin)
