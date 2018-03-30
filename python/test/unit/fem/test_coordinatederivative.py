"""Unit tests for coordinate derivative"""

# Copyright (C) 2011-2014 Garth N. Wells
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
# Modified by Florian Wechsung 2018

import pytest
import numpy as np
from dolfin import *


def test_first_shape_derivative():
    mesh = UnitSquareMesh(6, 6)
    n = FacetNormal(mesh)
    X = SpatialCoordinate(mesh)
    x, y = X
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = project(x*x+y*x+y*y+sin(x)+cos(x), V)
    dX = TestFunction(VectorFunctionSpace(mesh, "Lagrange", 1))

    J = u * u * dx
    computed = assemble(derivative(J, X, dX)).get_local()
    actual = assemble(u * u * div(dX) * dx).get_local()
    assert np.allclose(computed, actual, rtol=1e-14)    

    J = inner(grad(u), grad(u)) * dx
    computed = assemble(derivative(J, X)).get_local()
    dJdX = -2*inner(dot(grad(dX), grad(u)), grad(u)) * dx + inner(grad(u), grad(u)) * div(dX) * dx
    actual = assemble(dJdX).get_local()
    assert np.allclose(computed, actual, rtol=1e-14)    


    f = x * y * sin(x) * cos(y)
    J = f * dx
    computed = assemble(derivative(J, X, dX)).get_local()
    dJdX = div(f*dX) * dx
    actual = assemble(dJdX).get_local()
    assert np.allclose(computed, actual, rtol=1e-14)    

    J = f * ds
    computed = assemble(derivative(J, X, dX)).get_local()
    dJdX = inner(grad(f), dX) * ds + f * (div(dX) - inner(dot(grad(dX),n), n)) * ds
    actual = assemble(dJdX).get_local()
    assert np.allclose(computed, actual, rtol=1e-14)    

def test_mixed_derivatives():
    mesh = UnitSquareMesh(6, 6)
    X = SpatialCoordinate(mesh)
    x, y = X
    V = FunctionSpace(mesh, "CG", 1)
    u = project(x*x+y*x+y*y+sin(x)+cos(x), V)
    v = TrialFunction(V)
    dX = TestFunction(VectorFunctionSpace(mesh, "Lagrange", 1))
    dX_ = TrialFunction(VectorFunctionSpace(mesh, "Lagrange", 1))

    J = u * u * dx
    computed1 = assemble(derivative(derivative(J, X, dX), u)).array()
    computed2 = assemble(derivative(derivative(J, u), X, dX_)).array()
    actual = assemble(2 * u * v * div(dX) * dx).array()
    assert np.allclose(computed1, actual, rtol=1e-14)    
    assert np.allclose(computed2.T, actual, rtol=1e-14)    

    J = inner(grad(u), grad(u)) * dx
    computed1 = assemble(derivative(derivative(J, X, dX), u)).array()
    computed2 = assemble(derivative(derivative(J, u), X)).array()
    actual = assemble(2*inner(grad(u), grad(v)) * div(dX) * dx
                      - 2*inner(dot(grad(dX), grad(u)), grad(v)) * dx 
                      - 2*inner(grad(u), dot(grad(dX), grad(v))) * dx).array()
    assert np.allclose(computed1, actual, rtol=1e-14)    
    assert np.allclose(computed2.T, actual, rtol=1e-14)    

    
def test_second_shape_derivative():
    mesh = UnitSquareMesh(6, 6)
    V = FunctionSpace(mesh, "CG", 1)
    X = SpatialCoordinate(mesh)
    x, y = X
    u = project(x*x+y*x+y*y+sin(x)+cos(x), V)
    Z = VectorFunctionSpace(mesh, "Lagrange", 1)
    dX1 = TestFunction(Z)
    dX2 = TrialFunction(Z)

    J = u * u * dx
    computed = assemble(derivative(derivative(J, X), X)).array()
    actual = assemble(u * u * div(dX1) * div(dX2) * dx - u * u * tr(grad(dX1)*grad(dX2)) * dx).array()
    assert np.allclose(computed, actual, rtol=1e-14)    
    

if __name__ == "__main__":
    import os
    pytest.main(os.path.abspath(__file__))
