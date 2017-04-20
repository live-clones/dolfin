"""Unit tests for LocalPatchSolver"""

# Copyright (C) 2016-2017 Jan Blechta
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import pytest
import numpy

from dolfin import *
from dolfin_utils.test import pushpop_parameters, fixture


def test_interface(pushpop_parameters):
    # FIXME: Test for correct ghost mode in the solver class
    parameters["ghost_mode"] = "shared_vertex"
    mesh = UnitSquareMesh(32, 32)
    V = FunctionSpace(mesh, "P", 1)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = v*dx
    u = Function(V)
    x = u.vector()
    dofmap = V.dofmap()

    # Need a ghosted rhs
    # FIXME: Implement a check in LocalPatchSolver that rhs is ghosted,
    #        otherwise strange things happen
    b = Vector(mesh.mpi_comm())
    b_layout = b.factory().create_layout(b.rank())
    b_layout.init(b.mpi_comm(), [dofmap.index_map()], TensorLayout.Ghosts_GHOSTED)
    b.init(b_layout)
    assemble(L, tensor=b)

    with pytest.raises(RuntimeError):
        LocalPatchSolver([a, a])
    with pytest.raises(RuntimeError):
        LocalPatchSolver([])
    with pytest.raises(RuntimeError):
        LocalPatchSolver([a, a, a], [L, L])
    with pytest.raises(RuntimeError):
        LocalPatchSolver([a, a, a], [])

    LocalPatchSolver([a, a, a])
    LocalPatchSolver([a, a, a],
                     solver_type=LocalPatchSolver.SolverType_Cholesky,
                     bc_type=LocalPatchSolver.BCType_topological_zero)
    solver = LocalPatchSolver([a, a, a], [L, L, L],
                              solver_type=LocalPatchSolver.SolverType_Cholesky,
                              bc_type=LocalPatchSolver.BCType_topological_zero)

    solver.factorize()
    solver.clear_factorization()

    solver.solve_global_rhs(u)
    solver.solve_global_rhs([u, u, u])
    with pytest.raises(RuntimeError):
        solver.solve_global_rhs([u, u])
    with pytest.raises(RuntimeError):
        solver.solve_global_rhs([])

    solver.solve_local_rhs(u)
    solver.solve_local_rhs([u, u, u])
    with pytest.raises(RuntimeError):
        solver.solve_local_rhs([u, u])
    with pytest.raises(RuntimeError):
        solver.solve_local_rhs([])

    solver.solve_local(x, b, dofmap)
    solver.solve_local([x, x, x], [b, b, b], [dofmap, dofmap, dofmap])
    with pytest.raises(RuntimeError):
        solver.solve_local([x, x], [b, b, b], [dofmap, dofmap, dofmap])
    with pytest.raises(RuntimeError):
        solver.solve_local([], [b, b, b], [dofmap, dofmap, dofmap])
    with pytest.raises(RuntimeError):
        solver.solve_local([x, x, x], [b, b], [dofmap, dofmap, dofmap])
    with pytest.raises(RuntimeError):
        solver.solve_local([x, x, x], [b], [dofmap, dofmap, dofmap])
    with pytest.raises(RuntimeError):
        solver.solve_local([x, x, x], [b, b, b], [dofmap, dofmap])
    with pytest.raises(RuntimeError):
        solver.solve_local([x, x, x], [b, b, b], [])


@fixture
def mesh1d():
    n = 128
    parameters["ghost_mode"] = "shared_vertex"
    mesh = UnitIntervalMesh(n)
    parameters["ghost_mode"] = "none"
    return n, mesh

@fixture
def mesh2d():
    n = 32
    parameters["ghost_mode"] = "shared_vertex"
    mesh = UnitSquareMesh(n, n)
    parameters["ghost_mode"] = "none"
    return n, mesh

@fixture
def mesh3d():
    n = 8
    parameters["ghost_mode"] = "shared_vertex"
    mesh = UnitCubeMesh(n, n, n)
    parameters["ghost_mode"] = "none"
    return n, mesh

# pytest does not support fixtures in parametrize yet hence this metafixture
@pytest.fixture(params=['mesh1d', 'mesh2d', 'mesh3d'])
def mesh(request, pushpop_parameters):
        return request.getfuncargvalue(request.param)

def test_correctness(mesh):
    """Rudimentary test for correctness of LocalPatchSolver using simple
    examples with one dof per patch. Especially parallel correctness is
    to be excercised by this test."""
    n, mesh = mesh
    V = FunctionSpace(mesh, "P", 1)
    u, v = TrialFunction(V), TestFunction(V)
    a0 = u*v*dx
    a2 = inner(grad(u), grad(v))*dx
    L = v*dx
    u = Function(V)
    d = mesh.topology().dim()
    patches_per_cell = d + 1

    solver0 = LocalPatchSolver(patches_per_cell*[a0], patches_per_cell*[L],
                               solver_type=LocalPatchSolver.SolverType_Cholesky,
                               bc_type=LocalPatchSolver.BCType_topological_zero)
    solver2 = LocalPatchSolver(patches_per_cell*[a2], patches_per_cell*[L],
                               solver_type=LocalPatchSolver.SolverType_Cholesky,
                               bc_type=LocalPatchSolver.BCType_topological_zero)

    # Test solves with global rhs
    solver0.solve_global_rhs(u)
    assert numpy.allclose(u.vector().array(), 1/(d+1) / (2/(d+1) - 2/(d+2)))
    solver2.solve_global_rhs(u)
    assert numpy.allclose(u.vector().array(), 1/(d+1) / n**2 )

    # Need a ghosted rhs
    # FIXME: Implement a check in LocalPatchSolver that rhs is ghosted,
    #        otherwise strange things happen
    b = Vector(mesh.mpi_comm())
    #b_layout = b.factory().create_layout(b.rank())
    #b_layout.init(b.mpi_comm(), [V.dofmap().index_map()], TensorLayout.Ghosts_GHOSTED)
    #b.init(b_layout)
    assemble(L, tensor=b)

    # Test solves with assembled rhs
    solver0.solve_local(u.vector(), b, V.dofmap())
    assert numpy.allclose(u.vector().array(), 1/(d+1) / (2/(d+1) - 2/(d+2)))
    solver2.solve_local(u.vector(), b, V.dofmap())
    assert numpy.allclose(u.vector().array(), 1/(d+1) / n**2 )


def test_flux_reconstruction(pushpop_parameters):
    parameters["ghost_mode"] = "shared_vertex"
    mesh = UnitSquareMesh(32, 32)
    RT = FiniteElement('Raviart-Thomas', mesh.ufl_cell(), 1)
    DG = FiniteElement('Discontinuous Lagrange', mesh.ufl_cell(), 0)
    W = FunctionSpace(mesh, RT*DG)
    u, r = TrialFunctions(W)
    v, q = TestFunctions(W)
    a = ( inner(u, v) - inner(r, div(v)) - inner(q, div(u)) )*dx
    D = Expression(("x[1]", "-x[0]"), degree=1)
    L = ( inner(D, v) )*dx
    w = Function(W)

    solver = LocalPatchSolver([a, a, a], [L, L, L],
         solver_type=LocalPatchSolver.SolverType_LU,
         # FIXME: This is probably wrong!
         bc_type=LocalPatchSolver.BCType_topological_zero,
         # FIXME: This is probably wrong!
         constraint_type=LocalPatchSolver.ConstraintType_remove_last_dof)
    solver.solve_global_rhs(w)
