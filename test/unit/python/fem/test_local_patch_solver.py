#!/usr/bin/env py.test

"""Unit tests for LocalPatchSolver"""

# Copyright (C) 2016 Jan Blechta
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

import pytest
from dolfin import *


def test_interface():
    mesh = UnitSquareMesh(3, 3)
    V = FunctionSpace(mesh, "P", 1)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = v*dx
    u = Function(V)
    b = assemble(L)
    x = u.vector()
    dofmap = V.dofmap()

    with pytest.raises(RuntimeError):
        LocalPatchSolver([a, a])
    with pytest.raises(RuntimeError):
        LocalPatchSolver([])
    with pytest.raises(RuntimeError):
        LocalPatchSolver([a, a, a], [L, L])
    with pytest.raises(RuntimeError):
        LocalPatchSolver([a, a, a], [])

    LocalPatchSolver([a, a, a])
    solver = LocalPatchSolver([a, a, a], [L, L, L])

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
