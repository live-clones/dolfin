#!/usr/bin/env py.test

"""Unit tests for the FunctionSpaceProduct class"""

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
# Modified by Cecile Daversin-Catty 2017
#
# First added:  2017-07-25
# Last changed: 2017-07-25

import pytest
from dolfin import *
from ufl.log import UFLException

from dolfin_utils.test import fixture

@fixture
def mesh():
    return UnitCubeMesh(8, 8, 8)

@fixture
def V1(mesh):
    return FunctionSpace(mesh, 'CG', 1)

@fixture
def V2(mesh):
    return VectorFunctionSpace(mesh, 'CG', 1)

@fixture
def V3(mesh):
    return FunctionSpace(mesh, 'CG', 2)

@fixture
def W(V1,V2,V3):
    return FunctionSpaceProduct(V1,V2,V3)

@fixture
def W11(V1):
    return FunctionSpaceProduct(V1,V1)

@fixture
def W12(V1,V2):
    return FunctionSpaceProduct(V1,V2)

@fixture
def W13(V1,V3):
    return FunctionSpaceProduct(V1,V3)

@fixture
def f(W):
    return Function(W)

def test_python_interface(W, W11, W12, W13):
    assert isinstance(W11, FunctionSpaceProduct)
    assert isinstance(W12, FunctionSpaceProduct)
    assert isinstance(W13, FunctionSpaceProduct)
    assert isinstance(W, FunctionSpaceProduct)

    spaces_list = [W, W11, W12, W13]
    for space in spaces_list:
        assert len(space.sub_spaces()) == space.num_sub_spaces()
        for f in space.sub_spaces():
            assert isinstance(f, FunctionSpace)

def test_function(f):
    for i in f.num_sub_spaces():
        assert isinstance(f.sub(i).function_space(), FunctionSpace)
    
    
