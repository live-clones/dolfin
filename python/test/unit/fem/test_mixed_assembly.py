#!/usr/bin/env py.test

"""Unit tests for the extract_blocks function (based on FormSplitter ufl class)"""

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
# Modified by Cecile Daversin-Catty 2019
#
# First added:  2019-01-25
# Last changed: 2019-01-25

import pytest
from dolfin import *
from ufl.log import UFLException
from dolfin_utils.test import fixture

@fixture
def marker_2D2D():
    n = 20
    square = UnitSquareMesh(n, n)
    marker = MeshFunction("size_t", square, square.topology().dim(), 0)
    for c in cells(square):
        marker[c] = c.midpoint().x() < 0.5
    return marker

@fixture
def marker_3D2D():
    n = 10
    cube = UnitCubeMesh(n, n, n)
    marker = MeshFunction("size_t", cube, cube.topology().dim() - 1, 0)
    for f in facets(cube):
        marker[f] = 0.5 - DOLFIN_EPS < f.midpoint().z() < 0.5 + DOLFIN_EPS
    return marker

def mesh(marker, i):
    return MeshView.create(marker, i)

def space(mesh):
    return FunctionSpace(mesh, "Lagrange", 1)

def a(u,v):
    return inner(grad(u), grad(v))*dx

def L(f,v):
    return f*v*dx

def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

def boundary1(x):
    return x[0] < DOLFIN_EPS

def boundary2(x):
    return x[0] > 1.0 - DOLFIN_EPS

def params():
    return dict({"linear_solver":"direct"})

def test_mixed_assembly_2D2D(marker_2D2D, marker_3D2D):

    def _compare_solutions(marker):
        # Meshes
        _mesh1 = mesh(marker,0)
        _mesh2 = mesh(marker,1)
        # Spaces
        _space1 = space(_mesh1)
        _space2 = space(_mesh2)
        _space_mixed = FunctionSpaceProduct(_space1, _space2)
        # Trial functions
        _u1 = TrialFunction(_space1)
        _u2 = TrialFunction(_space2)
        (_u1_m, _u2_m) = TrialFunction(_space_mixed)
        # Test functions
        _v1 = TestFunction(_space1)
        _v2 = TestFunction(_space2)
        (_v1_m, _v2_m) = TestFunction(_space_mixed)
        # Bilinear forms
        _a1 = a(_u1,_v1)
        _a2 = a(_u2,_v2)
        _am = a(_u1_m, _v1_m) + a(_u2_m, _v2_m)
        # Linear forms
        f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
        _L1 = L(f,_v1)
        _L2 = L(f,_v2)
        _Lm = L(f,_v1_m) + L(f,_v2_m)
        # Solution
        sol1 = Function(_space1)
        sol2 = Function(_space2)
        sol = Function(_space_mixed)
        # BCs
        bc1 = DirichletBC(_space1, Constant(0.0), boundary1)
        bc2 = DirichletBC(_space2, Constant(0.0), boundary2)
        # Resolution
        solve(_a1 == _L1, sol1, bcs=bc1)
        solve(_a2 == _L2, sol2, bcs=bc2)
        solve(_am == _Lm, sol, bcs=[bc1,bc2], solver_parameters=params())
        
        assert len(sol1.vector()) == len(sol.sub(0).vector())
        for i in range(len(sol1.vector())):
            assert abs(sol1.vector()[i] - sol.sub(0).vector()[i]) < 1e-10
            
        assert len(sol2.vector()) == len(sol.sub(1).vector())
        for i in range(len(sol2.vector())):
            assert abs(sol2.vector()[i] - sol.sub(1).vector()[i]) < 1e-10

    _compare_solutions(marker_2D2D)
    # FIXME
    #_compare_solutions(marker_3D2D)
