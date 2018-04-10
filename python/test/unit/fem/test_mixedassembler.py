#!/usr/bin/env py.test

"""Unit tests for MixedAssembler"""

# Copyright (C) 2011 Johan Hake
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
#
# Modified by Cecile Daversin-Catty 2018
#
# First added:  2018-04-10
# Last changed: 2018-04-10

import pytest
from dolfin import *
from ufl.log import UFLException

from dolfin_utils.test import skip_in_parallel, fixture

@fixture
def cube():
    return UnitCubeMesh(2, 2, 2)

@fixture
def square():
    return UnitSquareMesh(5, 5)

@fixture
def exterior():
    return "exterior"

@fixture
def interior():
    return "interior"

@pytest.fixture(scope='module', params=range(2))
def meshes(cube, square, request):
    mesh = [cube, square]
    return mesh[request.param]

@pytest.fixture(scope='module', params=range(2))
def locations(interior, exterior, request):
    location = [interior, exterior]
    return location[request.param]
    
def test_offdiagblock_assembly(square, cube, interior, exterior):
    """Off-diagonal blocks assembly"""

    def submesh(mesh, location):
        marker = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
        if location == "interior":
            for f in facets(mesh):
                marker[f] = 0.5 - DOLFIN_EPS < f.midpoint().x() < 0.5 + DOLFIN_EPS
        elif location == "exterior":
            DomainBoundary().mark(marker, 1)
        return MeshView.create(marker, 1)

    def _check_scalar(mesh, location):
        V = FunctionSpace(mesh, 'CG', 2)
        Q = FunctionSpace(submesh(mesh, location), 'CG', 1)
        W = FunctionSpaceProduct(V, Q)

        f = Expression('x[0]+x[1]', degree=3)
        g = Expression('x[0]+3*x[1]', degree=3)

        (u, p) = TrialFunction(W)
        (v, q) = TestFunction(W)
        dx_ = Measure('dx', domain=W.sub_space(1).mesh())
        trace_form = inner(u, q)*dx_ + inner(u, v)*dx + inner(p, v)*dx_ + inner(p, q)*dx_
        rhs = inner(Constant(0), v)*dx + inner(Constant(0), q)*dx_

        AA, bb, _ = assemble_mixed_system(trace_form == rhs, Function(W))

        T = AA[2]
        Tt = AA[1]

        v = interpolate(f, V)
        q = interpolate(g, Q)

        Tv = Function(Q).vector()
        T.mult(v.vector(), Tv)
    
        qT = Function(V).vector()
        Tt.mult(q.vector(), qT)

        assert(abs(q.vector().inner(Tv) - v.vector().inner(qT)) <= 1e-12)
        
    def _check_vectorial(mesh, location):
        V = VectorFunctionSpace(mesh, 'CG', 2)
        Q = VectorFunctionSpace(submesh(mesh, location), 'CG', 1)
        W = FunctionSpaceProduct(V, Q)

        if mesh.topology().dim() == 2:
            f = Expression(('x[0]+x[1]', 'x[0]-x[1]'), degree=3)
            g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]'), degree=3)
        elif mesh.topology().dim() == 3:
            f = Expression(('x[0]+x[1]', 'x[0]-x[1]', 'x[0]+x[1]'), degree=3)
            g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]', 'x[0]+3*x[1]'), degree=3)

        (u, p) = TrialFunction(W)
        (v, q) = TestFunction(W)
        dx_ = Measure('dx', domain=W.sub_space(1).mesh())
        trace_form = inner(u, q)*dx_ + inner(u, v)*dx + inner(p, v)*dx_ + inner(p, q)*dx_

        if mesh.topology().dim() == 2:
            f = Expression(('x[0]+x[1]', 'x[0]-x[1]'), degree=3)
            g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]'), degree=3)
            rhs = inner(Constant((0, 0)), v)*dx + inner(Constant((0, 0)), q)*dx_
        elif mesh.topology().dim() == 3:
            f = Expression(('x[0]+x[1]', 'x[0]-x[1]', 'x[0]+x[1]'), degree=3)
            g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]', 'x[0]+3*x[1]'), degree=3)
            rhs = inner(Constant((0, 0, 0)), v)*dx + inner(Constant((0, 0, 0)), q)*dx_

        AA, bb, _ = assemble_mixed_system(trace_form == rhs, Function(W))

        T = AA[2]
        Tt = AA[1]

        v = interpolate(f, V)
        q = interpolate(g, Q)

        Tv = Function(Q).vector()
        T.mult(v.vector(), Tv)
    
        qT = Function(V).vector()
        Tt.mult(q.vector(), qT)

        assert(abs(q.vector().inner(Tv) - v.vector().inner(qT)) <= 1e-12)

    _check_scalar(square, exterior)
    _check_scalar(square, interior)
    _check_vectorial(square, exterior)
    _check_vectorial(square, interior)
    
    _check_scalar(cube, exterior)
    _check_scalar(cube, interior)
    _check_vectorial(cube, exterior)
    _check_vectorial(cube, interior)
