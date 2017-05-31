#!/usr/bin/env py.test

"""Unit tests for BoundingBoxTree with MeshEntity"""

# Copyright (C) 2013-2014 Anders Logg
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

from __future__ import print_function
import pytest
import numpy

from dolfin import Mesh, MeshEditor, MeshValueCollection, MeshFunction
from dolfin import BoundingBoxTree
from dolfin import UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh, UnitDiscMesh
from dolfin import Vertex, Facet, Cell
from dolfin import vertices, facets, cells
from dolfin import Point
from dolfin import MPI, mpi_comm_world
from dolfin_utils.test import skip_in_parallel

@skip_in_parallel
def test_compute_collisions_interval_interval_2d():
    meshes = (UnitSquareMesh(16, 16),
        UnitDiscMesh(mpi_comm_world(), 4, 1, 2))

    for mesh in meshes:
        # Build reference of facet neighbours
        reference = {}
        for f in facets(mesh):
            neighbours = set([])
            for v in vertices(f):
              neighbours |= set([neighbour_f.index()
                for neighbour_f in facets(v)])
            reference[f.index()] = neighbours

        bbox = BoundingBoxTree()
        bbox.build(mesh, 1)

        for f in facets(mesh):
            entities = bbox.compute_entity_collisions(f)
            assert set(entities) == reference[f.index()]

@skip_in_parallel
def test_compute_collisions_self_intersect_interval_interval_2d():

    # Facet index : intersecting facet indices
    reference = {0: set([0, 1, 3, 4, 5, 6, 7, 8]),
                 1: set([0, 1, 2, 5, 6, 7, 8]),
                 2: set([8, 1, 2, 9, 7]),
                 3: set([0, 10, 3, 4]),
                 4: set([0, 3, 4, 5, 6, 7, 10, 11, 12]),
                 5: set([0, 1, 4, 5, 6, 8, 10, 11, 12]),
                 6: set([0, 1, 4, 5, 6, 7, 8, 11, 12]),
                 7: set([0, 1, 2, 4, 6, 7, 8, 11, 12]),
                 8: set([0, 1, 2, 5, 6, 7, 8, 9, 12]),
                 9: set([8, 9, 2, 12]),
                 10: set([11, 10, 3, 4, 5]),
                 11: set([4, 5, 6, 7, 10, 11, 12]),
                 12: set([4, 5, 6, 7, 8, 9, 11, 12])}
    
    # 3 x 1 mesh folded back on itself for self intersection
    node_lx = 2.0
    node_rx = 1.0
    node_ry = 0.9
    editor = MeshEditor()
    mesh = Mesh()
    editor.open(mesh, 2, 2)
    editor.init_cells(6)
    editor.init_vertices(8)
    editor.add_cell(0, 0, 5, 4)
    editor.add_cell(1, 0, 1, 5)
    editor.add_cell(2, 1, 6, 5)
    editor.add_cell(3, 1, 2, 6)
    editor.add_cell(4, 2, 7, 6)
    editor.add_cell(5, 2, 3, 7)
    editor.add_vertex(0, 0.0, 0.0)
    editor.add_vertex(1, node_lx, 0.0)
    editor.add_vertex(2, node_rx, -1.0 + node_ry)
    editor.add_vertex(3, 3.0, 0.0)
    editor.add_vertex(4, 0.0, 1.0)
    editor.add_vertex(5, node_lx, 1.0)
    editor.add_vertex(6, node_rx, node_ry)
    editor.add_vertex(7, 3.0, 1.0)
    editor.close()
    mesh.init()

    # Set up mvc to specify indices of facets. We do this since the facet
    # indices are not guaranteed the same between DOLFIN versions.
    mvc = MeshValueCollection('size_t', mesh, 1)
    mvc.set_value(0, 0, 10)
    mvc.set_value(0, 1, 4)
    mvc.set_value(0, 2, 3)
    
    mvc.set_value(1, 0, 5)
    mvc.set_value(1, 1, 4)
    mvc.set_value(1, 2, 0)

    mvc.set_value(2, 0, 11)
    mvc.set_value(2, 1, 6)
    mvc.set_value(2, 2, 5)

    mvc.set_value(3, 0, 7)
    mvc.set_value(3, 1, 6)
    mvc.set_value(3, 2, 1)

    mvc.set_value(4, 0, 12)
    mvc.set_value(4, 1, 8)
    mvc.set_value(4, 2, 7)

    mvc.set_value(5, 0, 9)
    mvc.set_value(5, 1, 8)
    mvc.set_value(5, 2, 2)

    # Put the mvc into a FacetFunction which we can use as a lookup
    ff = MeshFunction('size_t', mesh, mvc)

    # Build bbox of facets
    bbox = BoundingBoxTree()
    bbox.build(mesh, 1)

    for f_idx in reference.keys():
        f = Facet(mesh, f_idx)
        entities = bbox.compute_entity_collisions(f)
        entities = [ff[Facet(mesh, idx)] for idx in entities]
        assert set(entities) == reference[ff[Facet(mesh, f_idx)]]