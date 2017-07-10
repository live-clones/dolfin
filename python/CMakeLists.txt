cmake_minimum_required(VERSION 3.5.0)

PROJECT(dolfin_testYYY)

set(CMAKE_BUILD_TYPE Release)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})
add_subdirectory(pybind11)


# FIXME: Put this into a FindFoo.cmake
# Get mpi4py include path
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import mpi4py, sys; sys.stdout.write(mpi4py.get_include())"
  OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR)

# Create the binding library
pybind11_add_module(cpp SHARED
  src/dolfin.cpp
  src/common/common.cpp
  src/experimental/experimental.cpp
  src/fem/fem.cpp
  src/function/function.cpp
  src/generation/generation.cpp
  src/geometry/geometry.cpp
  src/mesh/mesh.cpp
  src/io/io.cpp
  src/la/la.cpp
  src/parameter/parameter.cpp
  src/refinement/refinement.cpp
  )

target_include_directories(cpp PRIVATE ${MPI4PY_INCLUDE_DIR})

# Add DOLFIN libraries and other config
target_link_libraries(cpp PRIVATE pybind11::module dolfin)

#set_target_properties(dolfin_test PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
#  SUFFIX "${PYTHON_MODULE_EXTENSION}")
