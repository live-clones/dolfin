cmake_minimum_required(VERSION 3.5.0)

PROJECT(dolfin_pybind11)

# Configurable options for how we want to build
include(FeatureSummary)
option(DOLFIN_ENABLE_PETSC4PY "Compile with support for PETSc4py." ON)
option(DOLFIN_ENABLE_MPI4PY "Compile with support for MPI4py." ON)
add_feature_info(DOLFIN_ENABLE_PETSC4PY DOLFIN_ENABLE_PETSC4PY "Compile with support for PETSc4py.")
add_feature_info(DOLFIN_ENABLE_MPI4PY DOLFIN_ENABLE_MPI4PY "Compile with support for MPI4py.")

# Get options from env (difficult to pass them through setuptools)
if(DEFINED ENV{DOLFIN_ENABLE_PETSC4PY})
  set(DOLFIN_ENABLE_PETSC4PY $ENV{DOLFIN_ENABLE_PETSC4PY})
endif()
if(DEFINED ENV{DOLFIN_ENABLE_MPI4PY})
  set(DOLFIN_ENABLE_MPI4PY $ENV{DOLFIN_ENABLE_MPI4PY})
endif()

execute_process(
  COMMAND
    "${PYTHON_EXECUTABLE}" -c
    "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE _tmp_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")


# Find required packages
find_package(pybind11 REQUIRED CONFIG HINTS ${PYBIND11_DIR} ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Strict compiler flags
#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-Wall -Werror -pedantic" HAVE_PEDANTIC)
#if (HAVE_PEDANTIC)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pedantic")
#endif()

# Create the binding library
pybind11_add_module(cpp SHARED
  src/dolfin.cpp
  src/parameter.cpp
  src/adaptivity.cpp
  src/ale.cpp
  src/common.cpp
  src/fem.cpp
  src/function.cpp
  src/generation.cpp
  src/geometry.cpp
  src/graph.cpp
  src/log.cpp
  src/math.cpp
  src/mesh.cpp
  src/multistage.cpp
  src/ts.cpp
  src/io.cpp
  src/la.cpp
  src/nls.cpp
  src/refinement.cpp
  src/MPICommWrapper.cpp)

# Add DOLFIN libraries and other config
target_link_libraries(cpp PRIVATE pybind11::module dolfin)

# Add to CMake search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Check for petsc4py
set(PETSC4PY_FOUND FALSE)
if(DOLFIN_ENABLE_PETSC4PY)
  find_package(PETSc4py)
endif()
if (PETSC4PY_FOUND)
  target_include_directories(cpp PRIVATE ${PETSC4PY_INCLUDE_DIRS})
  target_compile_definitions(cpp PRIVATE HAS_PYBIND11_PETSC4PY)
endif()

# Check for mpi4py
set(MPI4PY_FOUND FALSE)
if(DOLFIN_ENABLE_MPI4PY)
  find_package(MPI4PY)
endif()
if (MPI4PY_FOUND)
  # Don't know, but the standard MPI4PY_INCLUDE_DIR doesn't seem to work?
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import mpi4py, sys; sys.stdout.write(mpi4py.get_include())"
    OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR1
  )
  target_include_directories(cpp PRIVATE ${MPI4PY_INCLUDE_DIR1})
  target_compile_definitions(cpp PRIVATE HAS_PYBIND11_MPI4PY)
endif()

# scikit-build specific steps
if (SKBUILD)
  # Add relative rpath so _gosscpp can find the Goss::goss libray
  # when the build is relocated
  if (APPLE)
    set_target_properties(cpp PROPERTIES INSTALL_RPATH "@loader_path/lib")
  else()
    set_target_properties(cpp PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
  endif()

  install(TARGETS cpp DESTINATION .)
endif()
