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
def one_element(): # Reference element
    mesh = UnitSquareMesh(1, 1, "left")
    cell_f = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    cell_f[0] = 1
    mesh = SubMesh(mesh, cell_f, 1)

    marker = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    for f in facets(mesh):
        marker[f] = 0.0 - DOLFIN_EPS < f.midpoint().x() < 0.0 + DOLFIN_EPS

    return (mesh, marker)

@fixture
def two_elements():
    mesh = UnitSquareMesh(1, 1, "left")

    marker = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    for f in facets(mesh):
        marker[f] = 0.5 - DOLFIN_EPS < f.midpoint().x() < 0.5 + DOLFIN_EPS and 0.5 - DOLFIN_EPS < f.midpoint().y() < 0.5 + DOLFIN_EPS

    return (mesh, marker)


## TODO : Add the case of an interface between two submeshes
## Same mesh as two_elements but each elements is a submesh
# @fixture
# def interface():
#     return "interface"

@pytest.fixture(scope='module', params=range(2))
def meshes(one_element, two_elements, request):
    mesh = [one_element, two_elements]
    return mesh[request.param]

def test_mixed_assembly(one_element, two_elements):
    """Off-diagonal blocks assembly"""

    ## Test for a range of various FE types
    def _check_scalar(case, order):
        mesh = case[0]
        submesh = MeshView.create(case[1], 1)

        V = FunctionSpace(mesh, 'CG', order[0])
        Q = FunctionSpace(submesh, 'CG', order[1])
        W = FunctionSpaceProduct(V, Q)

        f = Expression('x[0]+x[1]', degree=3)
        g = Expression('x[0]-x[1]', degree=3)

        (u, p) = TrialFunction(W)
        (v, q) = TestFunction(W)

        dx_ = Measure('dx', domain=W.sub_space(1).mesh())
        ## Reference value
        ref = assemble(inner(f, g)*dx_)

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

        vqT = v.vector().inner(qT)
        qTv = q.vector().inner(Tv)

        assert(abs(vqT - ref) <= 1e-12)
        assert(abs(qTv - ref) <= 1e-12)
        assert(abs(qTv - vqT) <= 1e-12)
        
    def _check_vectorial(case, order):
        mesh = case[0]
        submesh = MeshView.create(case[1], 1)

        V = VectorFunctionSpace(mesh, 'CG', order[0])
        Q = VectorFunctionSpace(submesh, 'CG', order[1])
        W = FunctionSpaceProduct(V, Q)

        f = Expression(('x[0]+x[1]', 'x[0]-x[1]'), degree=3)
        g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]'), degree=3)

        (u, p) = TrialFunction(W)
        (v, q) = TestFunction(W)
        dx_ = Measure('dx', domain=W.sub_space(1).mesh())
        ## Reference value
        ref = assemble(inner(f, g)*dx_)

        trace_form = inner(u, q)*dx_ + inner(u, v)*dx + inner(p, v)*dx_ + inner(p, q)*dx_

        f = Expression(('x[0]+x[1]', 'x[0]-x[1]'), degree=3)
        g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]'), degree=3)
        rhs = inner(Constant((0, 0)), v)*dx + inner(Constant((0, 0)), q)*dx_

        AA, bb, _ = assemble_mixed_system(trace_form == rhs, Function(W))

        T = AA[2]
        Tt = AA[1]

        v = interpolate(f, V)
        q = interpolate(g, Q)

        Tv = Function(Q).vector()
        T.mult(v.vector(), Tv)
        qT = Function(V).vector()
        Tt.mult(q.vector(), qT)

        vqT = v.vector().inner(qT)
        qTv = q.vector().inner(Tv)

        assert(abs(vqT - ref) <= 1e-12)
        assert(abs(qTv - ref) <= 1e-12)
        assert(abs(qTv - vqT) <= 1e-12)

    # Scalar case - One element + Lagrange mult on exterior boundary
    _check_scalar(one_element, (1,2)) # CG1 - CG2
    _check_scalar(one_element, (2,1)) # CG2 - CG1
    # Scalar case - Two elements + Lagrange mult on interior facet
    _check_scalar(two_elements, (1,2)) # CG1 - CG2
    _check_scalar(two_elements, (2,1)) # CG2 - CG1
    # Vectorial case - One element + Lagrange mult on exterior boundary
    _check_vectorial(one_element, (1,2)) # CG1 - CG2
    _check_vectorial(one_element, (2,1)) # CG2 - CG1
    # Vectorial case - Two elements + Lagrange mult on interior facet
    _check_vectorial(two_elements, (1,2)) # CG1 - CG2
    _check_vectorial(two_elements, (2,1)) # CG2 - CG1
