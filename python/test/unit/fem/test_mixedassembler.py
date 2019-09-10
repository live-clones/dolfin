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

@fixture
def two_elements_with_interface():
    mesh = UnitSquareMesh(1, 1)

    markerc = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    markerf = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    for c in cells(mesh):
        markerc[c] = c.midpoint().y() > c.midpoint().x()
    for f in facets(mesh):
        markerf[f] = abs(f.midpoint().x() - f.midpoint().y()) < 1e-10

    return (markerc, markerf)

@fixture
def unit_marker_2D2D():
    n = 20
    square = UnitSquareMesh(n, n)
    marker = MeshFunction("size_t", square, square.topology().dim(), 0)
    for c in cells(square):
        marker[c] = c.midpoint().x() < 0.5
    return marker

@fixture
def unit_marker_3D2D():
    n = 20
    cube = UnitCubeMesh(n, n, n)
    marker = MeshFunction("size_t", cube, cube.topology().dim() - 1, 0)
    for f in facets(cube):
        marker[f] = 0.5 - DOLFIN_EPS < f.midpoint().z() < 0.5 + DOLFIN_EPS
    return marker

def meshview(marker, i):
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
    return x[0] > 1.0 - DOLFIN_EPS

def boundary2(x):
    return x[0] < DOLFIN_EPS

def params():
    return dict({"linear_solver":"direct"})

@skip_in_parallel
def test_mixed_assembly(one_element, two_elements):
    """Off-diagonal blocks assembly"""

    ## Test for a range of various FE types
    def _check_scalar(case, order):
        mesh = case[0]
        submesh = MeshView.create(case[1], 1)

        V = FunctionSpace(mesh, 'CG', order[0])
        Q = FunctionSpace(submesh, 'CG', order[1])
        W = MixedFunctionSpace(V, Q)

        f = Expression('x[0]+x[1]', degree=3)
        g = Expression('x[0]-x[1]', degree=3)

        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)

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
        W = MixedFunctionSpace(V, Q)

        f = Expression(('x[0]+x[1]', 'x[0]-x[1]'), degree=3)
        g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]'), degree=3)

        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        dx_ = Measure('dx', domain=W.sub_space(1).mesh())
        ## Reference value
        ref = assemble(inner(f, g)*dx_)

        trace_form = inner(u, q)*dx_ + inner(u, v)*dx + inner(p, v)*dx_ + inner(p, q)*dx_
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

@skip_in_parallel
def test_mixed_assembly_interface(two_elements_with_interface):
    """Off-diagonal blocks assembly"""

    ## Test for a range of various FE types
    def _check_scalar(case, order):
        mesh0 = MeshView.create(case[0], 0)
        mesh1 = MeshView.create(case[0], 1)
        submesh = MeshView.create(case[1], 1)

        V0 = FunctionSpace(mesh0, 'CG', order[0])
        V1 = FunctionSpace(mesh1, 'CG', order[1])
        Q = FunctionSpace(submesh, 'CG', order[2])
        W = MixedFunctionSpace(V0, V1, Q)

        f = Expression('x[0]+x[1]', degree=3)
        g = Expression('x[0]-x[1]', degree=3)

        (u0, u1, p) = TrialFunctions(W)
        (v0, v1, q) = TestFunctions(W)

        dx_ = Measure('dx', domain=W.sub_space(2).mesh())
        # Reference value
        ref_fg = assemble(inner(f, g)*dx_)
        ref_ff = assemble(inner(f, f)*dx_)

        # Diagonal blocks
        trace_form = inner(u0, v0)*dx + inner(u1, v1)*dx + inner(p, q)*dx_
        # Mixed-dimensional coupling
        trace_form += inner(u0, q)*dx_ + inner(u1, q)*dx_ + inner(p, v0)*dx_ + inner(p, v1)*dx_
        # Integration of common
        trace_form += inner(u0, v1)*dx_ + inner(u1,v0)*dx_
        rhs = inner(Constant(0), v0)*dx + inner(Constant(0), v1)*dx + inner(Constant(0), q)*dx_

        AA, bb, _ = assemble_mixed_system(trace_form == rhs, Function(W))

        # Mixed-dimensional coupling
        T0 = AA[6] # u0*q
        T1 = AA[7] # u1*q
        Tt0 = AA[2] # p*v0
        Tt1 = AA[5] # p*v1
        # Integration of common interface
        I = AA[3] # u0*v1
        It = AA[1] # u1*v0

        v0 = interpolate(f, W.sub_space(0))
        v1 = interpolate(f, W.sub_space(1))
        q = interpolate(g, W.sub_space(2))

        Tv0 = Function(Q).vector()
        Tv1 = Function(Q).vector()
        T0.mult(v0.vector(), Tv0)
        T1.mult(v1.vector(), Tv1)
        qT0 = Function(V0).vector()
        qT1 = Function(V1).vector()
        Tt0.mult(q.vector(), qT0)
        Tt1.mult(q.vector(), qT1)

        Iv0 = Function(V1).vector()
        Iv1 = Function(V0).vector()
        I.mult(v0.vector(), Iv0)
        It.mult(v1.vector(), Iv1)

        v0qT0 = v0.vector().inner(qT0)
        v1qT1 = v1.vector().inner(qT1)
        qTv0 = q.vector().inner(Tv0)
        qTv1 = q.vector().inner(Tv1)

        v1Iv0 = v1.vector().inner(Iv0)
        v0Iv1 = v0.vector().inner(Iv1)

        assert(abs(v0qT0 - ref_fg) <= 1e-12)
        assert(abs(v1qT1 - ref_fg) <= 1e-12)
        assert(abs(qTv0 - ref_fg) <= 1e-12)
        assert(abs(qTv1 - ref_fg) <= 1e-12)
        assert(abs(qTv0 - v0qT0) <= 1e-12)
        assert(abs(qTv1 - v1qT1) <= 1e-12)

        assert(abs(v1Iv0 - ref_ff) <= 1e-12)
        assert(abs(v0Iv1 - ref_ff) <= 1e-12)

        # print("[scalar] v1Iv0 = ", v1Iv0, " while ref_ff = ", ref_ff)
        # print("[scalar] v0Iv1 = ", v0Iv1, " while ref_ff = ", ref_ff)


    def _check_vectorial(case, order):
        mesh0 = MeshView.create(case[0], 0)
        mesh1 = MeshView.create(case[0], 1)
        submesh = MeshView.create(case[1], 1)

        V0 = VectorFunctionSpace(mesh0, 'CG', order[0])
        V1 = VectorFunctionSpace(mesh1, 'CG', order[1])
        Q = VectorFunctionSpace(submesh, 'CG', order[2])
        W = MixedFunctionSpace(V0, V1, Q)

        f = Expression(('x[0]+x[1]', 'x[0]-x[1]'), degree=3)
        g = Expression(('x[0]+3*x[1]', 'x[0]-2*x[1]'), degree=3)

        (u0, u1, p) = TrialFunctions(W)
        (v0, v1, q) = TestFunctions(W)
        dx_ = Measure('dx', domain=W.sub_space(2).mesh())
        ## Reference value
        ref_fg = assemble(inner(f, g)*dx_)
        ref_ff = assemble(inner(f, f)*dx_)

        # Diagonal blocks
        trace_form = inner(u0, v0)*dx + inner(u1, v1)*dx + inner(p, q)*dx_
        # Mixed-dimensional coupling
        trace_form += inner(u0, q)*dx_ + inner(u1, q)*dx_ + inner(p, v0)*dx_ + inner(p, v1)*dx_
        # Integration of common interface
        trace_form += inner(u0, v1)*dx_ + inner(u1, v0)*dx_
        rhs = inner(Constant((0, 0)), v0)*dx + inner(Constant((0, 0)), v1)*dx + inner(Constant((0, 0)), q)*dx_

        AA, bb, _ = assemble_mixed_system(trace_form == rhs, Function(W))

        # Mixed-dimensional coupling
        T0 = AA[6]
        T1 = AA[7]
        Tt0 = AA[2]
        Tt1 = AA[5]
        # Integration of common interface
        I = AA[3] # u0*v1
        It = AA[1] # u1*v0

        v0 = interpolate(f, W.sub_space(0))
        v1 = interpolate(f, W.sub_space(1))
        q = interpolate(g, W.sub_space(2))

        Tv0 = Function(Q).vector()
        Tv1 = Function(Q).vector()
        T0.mult(v0.vector(), Tv0)
        T1.mult(v1.vector(), Tv1)
        qT0 = Function(V0).vector()
        qT1 = Function(V1).vector()
        Tt0.mult(q.vector(), qT0)
        Tt1.mult(q.vector(), qT1)

        Iv0 = Function(V1).vector()
        Iv1 = Function(V0).vector()
        I.mult(v0.vector(), Iv0)
        It.mult(v1.vector(), Iv1)

        v0qT0 = v0.vector().inner(qT0)
        v1qT1 = v1.vector().inner(qT1)
        qTv0 = q.vector().inner(Tv0)
        qTv1 = q.vector().inner(Tv1)

        v1Iv0 = v1.vector().inner(Iv0)
        v0Iv1 = v0.vector().inner(Iv1)

        assert(abs(v0qT0 - ref_fg) <= 1e-12)
        assert(abs(v1qT1 - ref_fg) <= 1e-12)
        assert(abs(qTv0 - ref_fg) <= 1e-12)
        assert(abs(qTv1 - ref_fg) <= 1e-12)
        assert(abs(qTv0 - v0qT0) <= 1e-12)
        assert(abs(qTv1 - v1qT1) <= 1e-12)

        assert(abs(v1Iv0 - ref_ff) <= 1e-12)
        assert(abs(v0Iv1 - ref_ff) <= 1e-12)

        # print("[vectorial] v1Iv0 = ", v1Iv0, " while ref_ff = ", ref_ff)
        # print("[vectorial] v0Iv1 = ", v0Iv1, " while ref_ff = ", ref_ff)

    # Scalar case - Two elements + Lagrange mult on interior facet
    _check_scalar(two_elements_with_interface, (1,1,1)) # CG1 - CG1 - CG2
    _check_scalar(two_elements_with_interface, (2,2,1)) # CG2 - CG2 - CG1
    _check_scalar(two_elements_with_interface, (1,1,2)) # CG1 - CG1- CG2
    _check_scalar(two_elements_with_interface, (2,2,2)) # CG2 - CG2- CG1
    _check_scalar(two_elements_with_interface, (1,1,2)) # CG2 - CG1- CG1
    _check_scalar(two_elements_with_interface, (1,2,1)) # CG1 - CG2- CG2
    _check_scalar(two_elements_with_interface, (2,1,1)) # CG1 - CG2- CG2
    # Vectorial case - Two elements + Lagrange mult on interior facet
    _check_vectorial(two_elements_with_interface, (1,1,2)) # CG1 - CG1 - CG2
    _check_vectorial(two_elements_with_interface, (2,2,1)) # CG2 - CG2 - CG1
    _check_vectorial(two_elements_with_interface, (1,1,2)) # CG1 - CG1- CG2
    _check_vectorial(two_elements_with_interface, (2,2,1)) # CG2 - CG2- CG1
    _check_vectorial(two_elements_with_interface, (2,1,1)) # CG2 - CG1- CG1
    _check_vectorial(two_elements_with_interface, (1,2,1)) # CG1 - CG2- CG2

@skip_in_parallel
def test_mixed_assembly_diag(unit_marker_2D2D, unit_marker_3D2D):
    def _compare_solutions(marker, boundaries):
        # Meshes
        _mesh1 = meshview(marker,0)
        _mesh2 = meshview(marker,1)
        # Spaces
        _space1 = space(_mesh1)
        _space2 = space(_mesh2)
        _space_mixed = MixedFunctionSpace(_space1, _space2)
        # Trial functions
        _u1 = TrialFunction(_space1)
        _u2 = TrialFunction(_space2)
        (_u1_m, _u2_m) = TrialFunctions(_space_mixed)
        # Test functions
        _v1 = TestFunction(_space1)
        _v2 = TestFunction(_space2)
        (_v1_m, _v2_m) = TestFunctions(_space_mixed)
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
        bc1 = DirichletBC(_space1, Constant(0.0), boundaries[0])
        bc2 = DirichletBC(_space2, Constant(0.0), boundaries[1])
        # Resolution
        solve(_a1 == _L1, sol1, bcs=bc1)
        solve(_a2 == _L2, sol2, bcs=bc2)
        solve(_am == _Lm, sol, bcs=[bc1, bc2], solver_parameters=params())

        sol1_m = sol.sub(0,deepcopy=True)
        sol2_m = sol.sub(1,deepcopy=True)
        assert len(sol1.vector()) == len(sol1_m.vector())
        for i in range(len(sol1.vector())):
            assert abs(sol1.vector()[i] - sol1_m.vector()[i]) < 1e-10

        assert len(sol2.vector()) == len(sol2_m.vector())
        for i in range(len(sol2.vector())):
            assert abs(sol2.vector()[i] - sol2_m.vector()[i]) < 1e-10


    _compare_solutions(unit_marker_2D2D, [boundary1, boundary2])
    _compare_solutions(unit_marker_3D2D, [boundary, boundary])
