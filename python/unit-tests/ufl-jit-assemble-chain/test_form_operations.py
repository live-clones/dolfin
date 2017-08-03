#!usr/bin/env py.test

"""Tests for DOLFIN integration of various form operations"""

# Copyright (C) 2011 Marie E. Rognes
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
import ufl
from dolfin import *

def test_lhs_rhs_simple():
    """Test taking lhs/rhs of DOLFIN specific forms (constants
    without cell). """

    mesh = RectangleMesh(Point(0, 0), Point(2, 1), 3, 5)
    V = FunctionSpace(mesh, "CG", 1)
    f = Constant(2.0)
    g = Constant(3.0)
    v = TestFunction(V)
    u = TrialFunction(V)

    F = inner(g*grad(f*v), grad(u))*dx + f*v*dx
    a, L = ufl.system(F)

    Fl = ufl.lhs(F)
    Fr = ufl.rhs(F)

    a0 = ufl.inner(grad(v), grad(u))*dx

    n = assemble(a).norm("frobenius")
    nl = assemble(Fl).norm("frobenius")
    n0 = 6.0*assemble(a0).norm("frobenius")

    assert round(n - n0, 7) == 0
    assert round(n - nl, 7) == 0
