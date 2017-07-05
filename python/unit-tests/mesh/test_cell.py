#!/usr/bin/env py.test

"""Unit tests for the Cell class"""

# Copyright (C) 2013 Anders Logg
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
# First added:  2013-04-18
# Last changed: 2014-05-30

import pytest
import numpy
#from dolfin import *
from dolfin_test.cpp.generation import (UnitIntervalMesh, UnitSquareMesh,
                                        UnitCubeMesh, UnitQuadMesh)
from dolfin_test.cpp.mesh import (Mesh, Cell, MeshEditor)
from dolfin_test.cpp.geometry import (Point)
import dolfin_test.cpp.MPI as MPI

#from dolfin_utils.test import skip_in_parallel, skip_in_release


#@skip_in_parallel
def test_distance_interval():

    mesh = UnitIntervalMesh(1)
    cell = Cell(mesh, 0)

    assert round(cell.distance(Point(-1.0)) - 1.0, 7) == 0
    assert round(cell.distance(Point(0.5)) - 0.0, 7) == 0


#@skip_in_parallel
def test_distance_triangle():

    mesh = UnitSquareMesh(1, 1)
    cell = Cell(mesh, 1)

    assert round(cell.distance(Point(-1.0, -1.0)) - numpy.sqrt(2), 7) == 0
    assert round(cell.distance(Point(-1.0, 0.5)) - 1, 7) == 0
    assert round(cell.distance(Point(0.5, 0.5)) - 0.0, 7) == 0


#@skip_in_parallel
def test_distance_tetrahedron():

    mesh = UnitCubeMesh(1, 1, 1)
    cell = Cell(mesh, 5)

    assert round(cell.distance(Point(-1.0, -1.0, -1.0))-numpy.sqrt(3), 7) == 0
    assert round(cell.distance(Point(-1.0, 0.5, 0.5)) - 1, 7) == 0
    assert round(cell.distance(Point(0.5, 0.5, 0.5)) - 0.0, 7) == 0


#@skip_in_release
#@skip_in_parallel
def test_issue_568():
    mesh = UnitSquareMesh(4, 4)
    cell = Cell(mesh, 0)

    # Should throw an error, not just segfault (only works in DEBUG mode!)
    with pytest.raises(RuntimeError):
        cell.facet_area(0)

    # Should work after initializing the connectivity
    mesh.init(2, 1)
    cell.facet_area(0)


def test_volume_quadrilateralR2():

    mesh = UnitQuadMesh(MPI.comm_self, 1, 1)
    cell = Cell(mesh, 0)

    assert cell.volume() == 1.0



def test_volume_quadrilateralR3_1():

    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, "quadrilateral", 2, 3)
    editor.init_vertices(4)
    editor.init_cells(1)
    editor.add_vertex(0, Point(0.0, 0.0, 0.0))
    editor.add_vertex(1, Point(1.0, 0.0, 0.0))
    editor.add_vertex(2, Point(0.0, 1.0, 0.0))
    editor.add_vertex(3, Point(1.0, 1.0, 0.0))
    editor.add_cell(0,numpy.array([0, 1, 2, 3],dtype=numpy.uintp))
    editor.close()
    mesh.init()
    cell = Cell(mesh, 0)

    assert cell.volume() == 1.0



def test_volume_quadrilateralR3_2():

    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, "quadrilateral", 2, 3)
    editor.init_vertices(4)
    editor.init_cells(1)
    editor.add_vertex(0, Point(0.0, 0.0, 0.0))
    editor.add_vertex(1, Point(0.0, 1.0, 0.0))
    editor.add_vertex(2, Point(0.0, 0.0, 1.0))
    editor.add_vertex(3, Point(0.0, 1.0, 1.0))
    editor.add_cell(0,numpy.array([0, 1, 2, 3],dtype=numpy.uintp))
    editor.close()
    mesh.init()
    cell = Cell(mesh, 0)

    assert cell.volume() == 1.0