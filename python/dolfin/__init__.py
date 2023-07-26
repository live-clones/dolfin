# -*- coding: utf-8 -*-
"""Main module for DOLFIN"""

# flake8: noqa

# Copyright (C) 2017 Chris N. Richardson and Garth N. Wells
#
# Distributed under the terms of the GNU Lesser Public License (LGPL),
# either version 3 of the License, or (at your option) any later
# version.

import sys

# Store dl open flags to restore them after import
stored_dlopen_flags = sys.getdlopenflags()

# Developer note: below is related to OpenMPI
# Fix dlopen flags (may need reorganising)
if "linux" in sys.platform:
    # FIXME: What with other platforms?
    try:
        from ctypes import RTLD_NOW, RTLD_GLOBAL
    except ImportError:
        RTLD_NOW = 2
        RTLD_GLOBAL = 256
    sys.setdlopenflags(RTLD_NOW | RTLD_GLOBAL)
del sys

# Reset dl open flags
# sys.setdlopenflags(stored_dlopen_flags)
# del sys

# Import cpp modules
from .cpp import __version__

from .cpp.common import (Variable, has_debug, has_hdf5, has_scotch,
                         has_hdf5_parallel, has_mpi, has_mpi4py,
                         has_petsc, has_petsc4py, has_parmetis, has_sundials,
                         has_slepc, has_slepc4py, git_commit_hash,
                         DOLFIN_EPS, DOLFIN_PI,  DOLFIN_EPS_LARGE,
                         TimingClear, TimingType,
                         timing, timings, list_timings, dump_timings_to_xml,
                         SubSystemsManager)

if has_hdf5():
    from .cpp.adaptivity import TimeSeries
    from .cpp.io import HDF5File

from .cpp.ale import ALE
from .cpp import MPI
from .cpp.function import (Expression, Constant, FunctionAXPY,
                           LagrangeInterpolator, FunctionAssigner, assign,
                           MultiMeshSubSpace)
from .cpp.fem import (FiniteElement, DofMap, Assembler, MultiMeshAssembler,
                      get_coordinates, create_mesh, set_coordinates,
                      vertex_to_dof_map, dof_to_vertex_map,
                      PointSource, DiscreteOperators,
                      LinearVariationalSolver,
                      NonlinearVariationalSolver,
                      MixedLinearVariationalSolver,
                      MixedNonlinearVariationalSolver,
                      SparsityPatternBuilder,
                      MultiMeshDirichletBC, adapt)

from .cpp.geometry import (BoundingBoxTree,
                           Point,
                           MeshPointIntersection,
                           intersect)
from .cpp.generation import (IntervalMesh, BoxMesh, RectangleMesh,
                             UnitDiscMesh,
                             UnitTriangleMesh, UnitCubeMesh,
                             UnitSquareMesh, UnitIntervalMesh,
                             SphericalShellMesh)
from .cpp.graph import GraphBuilder
from .cpp.io import File, XDMFFile, VTKFile
from .cpp.la import (list_linear_algebra_backends,
                     list_linear_solver_methods,
                     list_lu_solver_methods,
                     list_krylov_solver_methods,
                     list_krylov_solver_preconditioners,
                     has_linear_algebra_backend,
                     has_lu_solver_method,
                     has_krylov_solver_method,
                     has_krylov_solver_preconditioner,
                     linear_algebra_backends,
                     linear_solver_methods,
                     lu_solver_methods,
                     krylov_solver_methods,
                     krylov_solver_preconditioners,
                     normalize,
                     VectorSpaceBasis, in_nullspace,
                     residual)

if has_linear_algebra_backend('PETSc'):
    from .cpp.la import (PETScVector, PETScMatrix, PETScNestMatrix, PETScFactory,
                         PETScOptions, PETScLUSolver,
                         PETScKrylovSolver, PETScPreconditioner)
    from .cpp.fem import PETScDMCollection
    from .cpp.nls import (PETScSNESSolver, PETScTAOSolver, TAOLinearBoundSolver)

if has_linear_algebra_backend('Tpetra'):
    from .cpp.la import (TpetraVector, TpetraMatrix, TpetraFactory,
                         MueluPreconditioner, BelosKrylovSolver)

if has_slepc():
    from .cpp.la import SLEPcEigenSolver

from .cpp.la import (IndexMap, DefaultFactory, Matrix, Vector, Scalar,
                     EigenMatrix, EigenVector, EigenFactory, LUSolver,
                     KrylovSolver, TensorLayout, LinearOperator,
                     BlockMatrix, BlockVector)
from .cpp.la import GenericVector  # Remove when pybind11 transition complete
from .cpp.log import (info, Table, set_log_level, get_log_level, LogLevel,
                      Progress, begin, end, error, warning, set_log_active)
from .cpp.math import ipow, near, between
from .cpp.mesh import (Mesh, MeshTopology, MeshGeometry, MeshEntity,
                       MeshColoring, CellType, Cell, Facet, Face,
                       Edge, Vertex, cells, facets, faces, edges,
                       entities, vertices, SubDomain, BoundaryMesh,
                       MeshEditor, MeshQuality, SubMesh,
                       DomainBoundary, PeriodicBoundaryComputation,
                       MeshTransformation, SubsetIterator, MultiMesh, MeshView,
                       MeshPartitioning)

from .cpp.nls import (NonlinearProblem, NewtonSolver, OptimisationProblem)
from .cpp.refinement import refine, p_refine
from .cpp.parameter import Parameters, parameters
from .cpp.io import X3DOM, X3DOMParameters

if has_sundials():
    from .cpp.la import SUNDIALSNVector
    from .cpp.ts import CVode

# Import Python modules
from . import io
from . import la
from . import mesh
from . import parameter

from .common import timer
from .common.timer import Timer, timed
from .common.plotting import plot

from .fem.assembling import (assemble, assemble_system, assemble_multimesh, assemble_mixed,
                             SystemAssembler, assemble_local)
from .fem.form import Form
from .fem.norms import norm, errornorm
from .fem.dirichletbc import DirichletBC, AutoSubDomain
from .fem.multimeshdirichletbc import MultiMeshDirichletBC
from .fem.interpolation import interpolate
from .fem.projection import project
from .fem.solvers import LocalSolver
from .fem.solving import (solve, LinearVariationalProblem,
                          NonlinearVariationalProblem,
                          MixedLinearVariationalProblem,
                          MixedNonlinearVariationalProblem)
from .fem.solving import assemble_mixed_system, solve_mixed_system
from .fem.formmanipulations import (derivative, adjoint, increase_order, tear, extract_blocks)

# Need to be careful with other to avoid circular dependency
from .fem.adaptivesolving import (AdaptiveLinearVariationalSolver,
                                  AdaptiveNonlinearVariationalSolver)

from .function.multimeshfunctionspace import (MultiMeshFunctionSpace,
                                              MultiMeshVectorFunctionSpace,
                                              MultiMeshTensorFunctionSpace)
from .function.functionspace import (FunctionSpace,
                                     MixedFunctionSpace,
                                     VectorFunctionSpace, 
                                     TensorFunctionSpace)

from .function.function import Function
from .function.multimeshfunction import MultiMeshFunction
from .function.argument import (TestFunction, TrialFunction,
                                TestFunctions, TrialFunctions)
from .function.constant import Constant
from .function.specialfunctions import (MeshCoordinates, FacetArea, FacetNormal,
                                        CellVolume, SpatialCoordinate, CellNormal,
                                        CellDiameter, Circumradius,
                                        MinCellEdgeLength, MaxCellEdgeLength,
                                        MinFacetEdgeLength, MaxFacetEdgeLength)
from .function.expression import Expression, UserExpression, CompiledExpression

# experimental
from .jit.pybind11jit import compile_cpp_code

from .la import as_backend_type, la_index_dtype
from .mesh.ale import (compute_vertex_map, compute_edge_map,
                       init_parent_edge_indices)
from .mesh.meshfunction import (MeshFunction)
from .mesh.meshvaluecollection import MeshValueCollection
from .mesh.subdomain import CompiledSubDomain

from .multistage.multistagescheme import (RK4, CN2, CrankNicolson,
                                          ExplicitMidPoint,
                                          ESDIRK3, ESDIRK4,
                                          ForwardEuler, BackwardEuler)
from .multistage.multistagesolvers import PointIntegralSolver, RKSolver
from .multistage.rushlarsenschemes import RL1, RL2, GRL1, GRL2

# Import from ufl (2022.2 or earlier)
from ufl_legacy import (FiniteElement, TensorElement, VectorElement,
                 MixedElement, NodalEnrichedElement, rhs, lhs, conditional, le,
                 lt, ge, gt, split, cross, inner, dot, grad, nabla_grad, curl,
                 dx, div, Measure, det, pi, sin, cos, tan, acos, asin, atan,
                 ln, exp, sqrt, bessel_I, bessel_J, bessel_K,
                 bessel_Y, Dx, ds, dS, dP, dX, dC, dI, dO, interval, triangle,
                 tetrahedron, quadrilateral, hexahedron, avg, jump,
                 sym, tr, Identity, variable, diff, as_vector,
                 as_tensor, as_matrix, system, outer, dev, skew,
                 elem_mult, elem_div, elem_pow, elem_op, erf, inv)
from ufl_legacy.formoperators import action
