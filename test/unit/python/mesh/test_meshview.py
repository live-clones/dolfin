#!/usr/bin/env py.test

"""Unit tests for MeshViewMapping class"""

# Copyright (C) 2011 Garth N. Wells
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
# First added:  2011-02-26
# Last changed: 2014-05-30

import pytest
from dolfin import *
from dolfin_utils.test import fixture, skip_in_parallel


@fixture
def cube():
    return UnitCubeMesh(5, 5, 5)


@fixture
def square():
    return UnitSquareMesh(5, 5)


@pytest.fixture(scope='module', params=range(2))
def meshes(cube, square, request):
    mesh = [cube, square]
    return mesh[request.param]

def test_make_facets_view(square, cube):
    """Create view of facets"""

    def _check_facets_view(mesh):
        marker = FacetFunction("size_t", mesh, 0)
        D = mesh.topology().dim()
        mesh.init(D-1, D)

        c = 0
        for f in facets(mesh):
            if f.num_global_entities(D) == 1:
                marker[f] = 1
                c += 1

        m2 = MeshViewMapping.create_from_marker(marker, 1)

        assert m2.num_cells() == c
        assert m2.num_cells() == len(m2.topology().mapping.cell_map())
        assert m2.num_vertices() == len(m2.topology().mapping.vertex_map())

    _check_facets_view(square)
    _check_facets_view(cube)


def test_make_cells_view(square, cube):
    """Create view of cells"""

    def _check_cells_view(mesh):
        marker = CellFunction("size_t", mesh, 0)

        ct = 0
        for c in cells(mesh):
            if c.midpoint().x() < 0.5:
                marker[c] = 1
                ct += 1

        m2 = MeshViewMapping.create_from_marker(marker, 1)

        assert m2.num_cells() == ct
        assert m2.num_cells() == len(m2.topology().mapping.cell_map())
        assert m2.num_vertices() == len(m2.topology().mapping.vertex_map())

    _check_cells_view(square)
    _check_cells_view(cube)


def test_make_edges_view(cube):
    """Create view of edges"""

    marker = EdgeFunction("size_t", cube, 0)
    cube.init(1, 3)

    c = 0
    for e in edges(cube):
        if (e.midpoint().x() == 0 or e.midpoint().x() == 1):
            marker[e] = 1
            c += 1
        elif (e.midpoint().y() == 0 or e.midpoint().y() == 1):
            marker[e] = 1
            c += 1
        elif (e.midpoint().z() == 0 or e.midpoint().z() == 1):
            marker[e] = 1
            c += 1

    m2 = MeshViewMapping.create_from_marker(marker, 1)

    assert m2.num_cells() == c
    assert m2.num_cells() == len(m2.topology().mapping.cell_map())
    assert m2.num_vertices() == len(m2.topology().mapping.vertex_map())
