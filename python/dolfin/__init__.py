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
from .cpp.common import Variable, has_debug, has_hdf5, has_hdf5_parallel, \
    has_mpi, has_petsc, has_slepc, git_commit_hash
from .cpp import MPI
from .cpp.function import Expression, Constant, Function, interpolate
from .cpp.fem import FiniteElement, DofMap, Assembler
from .cpp.geometry import BoundingBoxTree, Point, MeshPointIntersection, intersect
from .cpp.generation import IntervalMesh, UnitIntervalMesh, \
    UnitSquareMesh, UnitCubeMesh, BoxMesh, \
    RectangleMesh, UnitQuadMesh
from .cpp.io import XDMFFile
from .cpp.la import has_linear_algebra_backend
if has_linear_algebra_backend('PETSc'):
    from .cpp.la import PETScVector, PETScMatrix
from .cpp.la import Matrix, Vector, EigenMatrix, EigenVector, LUSolver, KrylovSolver
from .cpp.mesh import Mesh, MeshTopology, MeshGeometry, MeshEntity, Cell, Facet, Face, Edge, Vertex, \
    cells, facets, faces, edges, vertices, SubDomain, BoundaryMesh
from .cpp.parameter import parameters
from .cpp.refinement import refine

# python modules
from .fem.form import Form
from .fem.dirichletbc import DirichletBC, CompiledSubDomain
from .function.functionspace import FunctionSpace
from .function.constant import Constant
from .function.expression import CompiledExpression, UserExpression
from .mesh.meshfunction import MeshFunction
