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
# Modified by Cecile Daversin-Catty 2017
#
# First added:  2017-07-25
# Last changed: 2017-07-25

import pytest
from dolfin import *
from ufl.log import UFLException

from dolfin_utils.test import fixture

from dolfin.fem.assembling import _create_dolfin_form

@fixture
def marker():
    cube = UnitCubeMesh(8, 8, 8)
    marker = CellFunction("size_t", cube, 0)
    for c in cells(cube):
        marker[c] = c.midpoint().x() < 0.5
    return marker

@fixture
def mesh1(marker):
    return MeshViewMapping.create_from_marker(marker, 1)

@fixture
def mesh2(marker):
    return MeshViewMapping.create_from_marker(marker, 0)

@fixture
def space(mesh1,mesh2):
    W1 = FunctionSpace(mesh1, "Lagrange", 1)
    W2 = FunctionSpace(mesh2, "Lagrange", 2)
    return FunctionSpaceProduct( W1, W2 )

@fixture
def u(space):
    return TrialFunction(space)

@fixture
def v(space):
    return TestFunction(space)

@fixture
def a(u,v):
    return inner(grad(u[0]), grad(v[0]))*dx + inner(grad(u[1]), grad(v[1]))*dx

# TODO : bilinear form with interface terms (non-diagonal blocks)
# @fixture
# def a2(u,v):
#     return ...

@fixture
def f():
    return Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)

@fixture
def L(f,v):
    return f*v[0]*dx + f*v[1]*dx

@fixture
def a_blocks(a):
    return extract_blocks(a)

@fixture
def L_blocks(L):
    return extract_blocks(L)

def test_size(space, a_blocks, L_blocks):
    assert len(L_blocks) == space.num_sub_spaces()
    assert len(a_blocks) == len(L_blocks)*len(L_blocks)

def test_indices(space, a, L, a_blocks, L_blocks):
    for i in range(space.num_sub_spaces()):   
        assert extract_blocks(L,i) == L_blocks[i]
        for j in range(space.num_sub_spaces()):
            assert extract_blocks(a,i,j) == a_blocks[i*space.num_sub_spaces()+j]

def test_spaces(space, a_blocks, L_blocks):
    for i in range(space.num_sub_spaces()):
        L_form = _create_dolfin_form(L_blocks[i])
        for j in range(space.num_sub_spaces()):
            if i == j : ## TEMPORARY : Do not conider non-diagonal blocks for now
                a_form = _create_dolfin_form(a_blocks[i*space.num_sub_spaces() + j])
                assert a_form.function_space(0) == L_form.function_space(0)
            # TODO : When non-diagonal blocks will be ok
            # if i < space.num_sub_spaces()-1 :
            #     a_form_next_row = _create_dolfin_form(a_blocks[(i+1)*space.num_sub_spaces() + j])
            #     assert a_form_next_row.function_space(1) == a_form.function_space(1)

        
        
        
