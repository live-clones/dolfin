# -*- coding: utf-8 -*-
"""Main module for DOLFIN"""

import sys

# Store dl open flags to restore them after import
stored_dlopen_flags = sys.getdlopenflags()

# Fix dlopen flags (may need reorganising)
import sys
if "linux" in sys.platform:
    # FIXME: What with other platforms?
    try:
        from ctypes import RTLD_NOW, RTLD_GLOBAL
    except ImportError:
        RTLD_NOW = 2
        RTLD_GLOBAL = 256
    sys.setdlopenflags(RTLD_NOW | RTLD_GLOBAL)
del sys

#import dolfin.cpp

# Reset dl open flags
#import sys
#sys.setdlopenflags(stored_dlopen_flags)
#del sys

# cpp modules
from .cpp.common import (Variable, has_debug, has_hdf5,
                         has_hdf5_parallel, has_mpi, has_petsc,
                         has_slepc, git_commit_hash, DOLFIN_EPS,
                         DOLFIN_PI, TimingClear, TimingType, timing)

if has_hdf5():
    from .cpp.adaptivity import TimeSeries
    from .cpp.io import HDF5File

from .cpp.ale import ALE
from .cpp import MPI
from .cpp.function import Expression, Constant #, interpolate
from .cpp.fem import (FiniteElement, DofMap, Assembler, SystemAssembler, get_coordinates,
                      set_coordinates, vertex_to_dof_map, dof_to_vertex_map, PointSource,
                      DiscreteOperators, assemble_local)
from .cpp.geometry import BoundingBoxTree, Point, MeshPointIntersection, intersect
from .cpp.generation import (IntervalMesh, BoxMesh, RectangleMesh, UnitDiscMesh, UnitQuadMesh, UnitHexMesh,
                             UnitCubeMesh, UnitSquareMesh, UnitIntervalMesh)
from .cpp.graph import GraphBuilder
from .cpp.io import File, XDMFFile, VTKFile
from .cpp.la import (has_linear_algebra_backend,
                     linear_algebra_backends, has_krylov_solver_method,
                     has_krylov_solver_preconditioner)

if has_linear_algebra_backend('PETSc'):
    from .cpp.la import PETScVector, PETScMatrix, PETScFactory
    from .cpp.fem import PETScDMCollection
    from .cpp.nls import PETScSNESSolver, PETScTAOSolver, TAOLinearBoundSolver

from .cpp.la import (IndexMap, DefaultFactory, Matrix, Vector, EigenMatrix,
                     EigenVector, EigenFactory, LUSolver, KrylovSolver)
from .cpp.log import info
from .cpp.math import ipow, near, between
from .cpp.mesh import (Mesh, MeshTopology, MeshGeometry, MeshEntity,
                       Cell, Facet, Face, Edge, Vertex, cells,
                       facets, faces, edges, entities,
                       vertices, SubDomain, BoundaryMesh,
                       MeshEditor, MultiMesh, MeshQuality,
                       SubMesh)
from .cpp.nls import NonlinearProblem, NewtonSolver
from .cpp.refinement import refine

from .cpp.parameter import Parameters, parameters

from .cpp.io import X3DOM, X3DOMParameters

# Python modules
from . import mesh
from . import parameter

from .common import timer
from .common.timer import Timer, timed
from .fem.assembling import assemble, assemble_system
from .fem.form import Form
from .fem.dirichletbc import DirichletBC, CompiledSubDomain
from .fem.interpolation import interpolate
from .function.functionspace import FunctionSpace, VectorFunctionSpace #, TensorFunctionSpace
from .function.function import Function
from .function.argument import TestFunction, TrialFunction, TestFunctions, TrialFunctions
from .function.constant import Constant
from .function.specialfunctions import FacetNormal, CellSize, SpatialCoordinate
from .function.expression import CompiledExpression, UserExpression
# from .mesh.mesh import UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh
from .la import as_backend_type
from .mesh.meshfunction import (MeshFunction, CellFunction,
                                FacetFunction, FaceFunction, EdgeFunction, VertexFunction)
from .mesh.meshvaluecollection import MeshValueCollection


# ufl
from ufl import (FiniteElement, VectorElement, MixedElement,
                 inner, dot, grad, dx, div,
                 ds, dS, triangle, tetrahedron, avg, jump)
from ufl.formoperators import action

# FIXME
def has_petsc4py():
    return False


def mpi_comm_world():
    return MPI.comm_world
