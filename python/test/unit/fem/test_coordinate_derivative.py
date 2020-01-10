"""Unit tests for coordinate derivative"""

# Copyright (C) 2018 Florian Wechsung
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
# Modified by Florian Wechsung 2018 and Jørgen S. Dokken 2020

import pytest
import numpy as np
from dolfin import *
from ufl import replace
from ufl.log import UFLException
from dolfin_utils.test import skip_in_parallel

def test_first_shape_derivative():
    mesh = UnitSquareMesh(6, 6)
    n = FacetNormal(mesh)
    X = SpatialCoordinate(mesh)
    x, y = X
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = project(x*x+y*x+y*y+sin(x)+cos(x), V)
    dX = TestFunction(VectorFunctionSpace(mesh, "Lagrange", 1))

    def test_first(J, dJ):
        computed = assemble(derivative(J, X)).get_local()
        actual = assemble(dJ).get_local()
        assert np.allclose(computed, actual, rtol=1e-14)

    Ja = u * u * dx
    dJa = u * u * div(dX) * dx
    test_first(Ja, dJa)

    Jb = inner(grad(u), grad(u)) * dx
    dJb = -2*inner(dot(grad(dX), grad(u)), grad(u)) * dx + inner(grad(u), grad(u)) * div(dX) * dx
    test_first(Jb, dJb)

    f = x * y * sin(x) * cos(y)
    Jc = f * dx
    dJc = div(f * dX) * dx
    test_first(Jc, dJc)

    Jd = f * ds
    dJd = inner(grad(f), dX) * ds \
        + f * (div(dX) - inner(dot(grad(dX), n), n)) * ds
    test_first(Jd, dJd)

    J = Ja + Jb + Jc + Jd
    dJ = dJa + dJb + dJc + dJd
    test_first(J, dJ)


def test_mixed_derivatives():
    mesh = UnitSquareMesh(6, 6)
    X = SpatialCoordinate(mesh)
    x, y = X
    V = FunctionSpace(mesh, "CG", 1)
    u = project(x * x + y * x + y * y + sin(x) + cos(x), V)
    v = TrialFunction(V)
    dX = TestFunction(VectorFunctionSpace(mesh, "Lagrange", 1))
    dX_ = TrialFunction(VectorFunctionSpace(mesh, "Lagrange", 1))


    def test_mixed(J, dJ_manual):
        computed1 = assemble(derivative(derivative(J, X, dX), u)).norm("frobenius")
        computed2 = assemble(derivative(derivative(J, u), X, dX_)).norm("frobenius")
        computed3 = assemble(derivative(derivative(J, X), u)).norm("frobenius")
        computed4 = assemble(derivative(derivative(J, u), X)).norm("frobenius")
        actuala = assemble(dJ_manual).norm("frobenius")
        assert np.isclose(computed1, actuala, rtol=1e-14)
        assert np.isclose(computed2, actuala, rtol=1e-14)
        assert np.isclose(computed3, actuala, rtol=1e-14)
        assert np.isclose(computed4, actuala, rtol=1e-14)

    Ja = u * u * dx
    dJa = 2 * u * v * div(dX) * dx
    test_mixed(Ja, dJa)

    Jb = inner(grad(u), grad(u)) * dx
    dJb = 2*inner(grad(u), grad(v)) * div(dX) * dx \
        - 2*inner(dot(nabla_grad(dX), grad(u)), grad(v)) * dx \
        - 2*inner(grad(u), dot(nabla_grad(dX), grad(v))) * dx
    test_mixed(Jb, dJb)

    J = Ja + Jb
    dJ = dJa + dJb
    test_mixed(J, dJ)



def test_second_shape_derivative():
    mesh = UnitSquareMesh(6, 6)
    V = FunctionSpace(mesh, "CG", 1)
    X = SpatialCoordinate(mesh)
    x, y = X
    u = project(x * x + y * x + y * y + sin(x) + cos(x), V)
    Z = VectorFunctionSpace(mesh, "Lagrange", 1)
    dX1 = TestFunction(Z)
    dX2 = TrialFunction(Z)

    def test_second(J, ddJ):
        computed = assemble(derivative(derivative(J, X, dX1), X, dX2)).norm("frobenius")
        actual = assemble(ddJ).norm("frobenius")
        assert np.isclose(computed, actual, rtol=1e-14)

    Ja = u * u * dx
    ddJa = u * u * div(dX1) * div(dX2) * dx - u * u * tr(grad(dX1)*grad(dX2)) * dx
    test_second(Ja, ddJa)

    Jb = inner(grad(u), grad(u)) * dx
    ddJb = 2*inner(dot(dot(nabla_grad(dX2), nabla_grad(dX1)), grad(u)), grad(u)) * dx \
        + 2*inner(dot(nabla_grad(dX1), dot(nabla_grad(dX2), grad(u))), grad(u)) * dx \
        + 2*inner(dot(nabla_grad(dX1), grad(u)), dot(nabla_grad(dX2), grad(u))) * dx \
        - 2*inner(dot(nabla_grad(dX2), grad(u)), grad(u)) * div(dX1) * dx \
        - inner(grad(u), grad(u)) * tr(nabla_grad(dX1)*nabla_grad(dX2)) * dx \
        - 2*inner(dot(nabla_grad(dX1), grad(u)), grad(u)) * div(dX2) * dx \
        + inner(grad(u), grad(u)) * div(dX1) * div(dX2) * dx
    test_second(Jb, ddJb)

    test_second(Ja+Jb, ddJa + ddJb)


# In parallel a RuntimeError is thrown instead of the specific
# UFLException, so we skip that case.
@skip_in_parallel
def test_integral_scaling_edge_case():
    mesh = UnitSquareMesh(6, 6)
    X = SpatialCoordinate(mesh)
    V = FunctionSpace(mesh, "CG", 1)
    u = Function(V)

    J = u * u * dx

    with pytest.raises(UFLException):
        assemble(Constant(2.0) * derivative(J, X))
    with pytest.raises(UFLException):
        assemble(derivative(Constant(2.0) * derivative(J, X), X))
    with pytest.raises(UFLException):
        assemble(Constant(2.0) * derivative(derivative(J, X), X))


def test_different_quadratures():
    """
    Checking that expressions with different quadrature rules
    can be assembled as a sum, using the higher order quadrature for all
    terms
    """
    mesh = UnitIntervalMesh(10)
    x = SpatialCoordinate(mesh)
    f = x[0]**2
    V = FunctionSpace(mesh, "CG", 1)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v)) * dx + inner(u, v) * dx
    l = f * v * dx
    bc = DirichletBC(V, Constant(1), "on_boundary")
    uh = Function(V)
    solve(a==l, uh, bcs=bc)
    F = a - l
    J = uh**2 * dx

    # Adjoint equation disregarding Dirichlet BC adjoint
    dJdu = derivative(J, uh)
    dFdu = derivative(action(F, uh), uh)
    bc_hom = DirichletBC(V, Constant(0), "on_boundary")
    lmbd = Function(V)
    solve(dFdu==-dJdu, lmbd, bc_hom)
    F = replace(F, {u: uh})
    dJdu = replace(dJdu, {u: uh})

    # Parts of second order adjoint equation:
    S = VectorFunctionSpace(mesh, "CG", 1)
    s = Function(S)
    F = replace(F, {v: lmbd})
    ds = TestFunction(S)
    dFdm = derivative(F, x, ds)
    d2Fdm2 = derivative(dFdm, x, s)
    soas_separated = assemble(d2Fdm2) + assemble(dFdm)

    soa = d2Fdm2 + dFdm
    soas = assemble(soa)
    assert np.isclose(soas_separated.norm("l2"), soas.norm("l2"))

    # 2D test for triangles with quadrature degree 3 and 4
    mesh = UnitSquareMesh(5, 5)
    x = SpatialCoordinate(mesh)
    V = FunctionSpace(mesh, "CG", 1)
    u0 = cos(pi*x[0])*sin(pi*x[1])

    u, v = TrialFunction(V), TestFunction(V)
    F = inner(u - u0, v)*dx
    uh = Function(V)
    solve(lhs(F) == rhs(F), uh)
    J = uh**2*dx

    # adjoint eq disregarding Dirichlet BC
    dJdu = derivative(J, uh)
    dFdu = derivative(action(F, uh), uh)
    bc_hom = DirichletBC(V, Constant(0), "on_boundary")
    lmbd = Function(V)
    solve(dFdu==-dJdu, lmbd, bc_hom)
    F = replace(F, {u:uh})
    dJdu = replace(dJdu, {u: uh})
    # Parts of second order adjoint eq:
    S = VectorFunctionSpace(mesh, "CG", 1)
    s = Function(S)
    F = replace(F, {v: lmbd})
    ds = TestFunction(S)
    dFdm = derivative(F, x, ds)
    d2Fdm2 = derivative(dFdm, x, s)
    d2Fdudm = derivative(dFdm, uh, lmbd)
    soas_separated = assemble(d2Fdm2) + assemble(dFdm) + assemble(d2Fdudm)
    soa = d2Fdm2 + d2Fdudm + dFdm
    soas = assemble(soa)
    assert np.isclose(soas_separated.norm("l2"), soas.norm("l2"))


def test_unique_tables():
    """
    Checking that the unique table names in FFC are cached correctly
    """
    mesh = UnitIntervalMesh(10)
    x = SpatialCoordinate(mesh)
    f = x[0]**2
    V = FunctionSpace(mesh, "CG", 1)
    f = project(f, V)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v)) * dx + inner(u, v) * dx
    l = f * v * dx
    bc = DirichletBC(V, Constant(1), "on_boundary")
    uh = Function(V)
    solve(a == l, uh, bcs=bc)
    F = a - l
    J = uh**2 * dx

    # Adjoint eq disregarding Dirichlet BC
    dJdu = derivative(J, uh)
    dFdu = derivative(action(F, uh), uh)
    bc_hom = DirichletBC(V, Constant(0), "on_boundary")
    lmbd = Function(V)
    solve(dFdu == -dJdu, lmbd, bc_hom)
    F = replace(F, {u: uh})
    dJdu = replace(dJdu, {u: uh})

    # Parts of second order adjoint eq:
    S = VectorFunctionSpace(mesh, "CG", 1)
    s = Function(S)
    F = replace(F, {v: lmbd})
    ds = TestFunction(S)
    dFdm = derivative(F, x, ds)
    d2Fdm2 = derivative(dFdm, x, s)
    soas_separated = assemble(d2Fdm2) + assemble(dFdm)

    soa = d2Fdm2 + dFdm
    soas = assemble(soa)
    assert np.isclose(soas_separated.norm("l2"), soas.norm("l2"))


if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
