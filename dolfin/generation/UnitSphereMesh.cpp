// Copyright (C) 2017 Chris Richardson
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//

#include <dolfin/common/MPI.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/refinement/refine.h>
#include "UnitSphereMesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void UnitSphereMesh::build(Mesh& mesh, std::size_t nrefine, std::size_t degree)
{
  MeshEditor editor;
  const std::size_t tdim = 3;
  const std::size_t gdim = 3;

  dolfin_assert(degree > 0 and degree < 3);

  Mesh base_mesh(mesh.mpi_comm());
  editor.open(base_mesh, tdim, gdim, 1);

  editor.init_vertices_global(13, 13);

  const double l0 = 2.0/(sqrt(10.0 + 2.0*sqrt(5.0)));
  const double l1 = l0*(1.0 + sqrt(5.0))/2.0;

  // Generate an icosahedron

  editor.add_vertex(0,  Point(  0,  l0, l1));
  editor.add_vertex(1,  Point(  0,  l0, -l1));
  editor.add_vertex(2,  Point(  0, -l0, -l1));
  editor.add_vertex(3,  Point(  0, -l0, l1));
  editor.add_vertex(4,  Point( l1,   0, l0));
  editor.add_vertex(5,  Point(-l1,   0, l0));
  editor.add_vertex(6,  Point(-l1,   0, -l0));
  editor.add_vertex(7,  Point( l1,   0, -l0));
  editor.add_vertex(8,  Point( l0,  l1, 0));
  editor.add_vertex(9,  Point( l0, -l1, 0));
  editor.add_vertex(10, Point(-l0, -l1, 0));
  editor.add_vertex(11, Point(-l0,  l1, 0));
  editor.add_vertex(12, Point(  0,   0, 0));

  editor.init_cells_global(20, 20);

  editor.add_cell(0, 0, 4, 8, 12);
  editor.add_cell(1, 0, 5, 11, 12);
  editor.add_cell(2, 1, 6, 11, 12);
  editor.add_cell(3, 1, 7, 8, 12);
  editor.add_cell(4, 2, 6, 10, 12);
  editor.add_cell(5, 2, 7, 9, 12);
  editor.add_cell(6, 3, 4, 9, 12);
  editor.add_cell(7, 3, 5, 10, 12);
  editor.add_cell( 8, 0, 3, 4, 12);
  editor.add_cell( 9, 0, 3, 5, 12);
  editor.add_cell(10, 1, 2, 6, 12);
  editor.add_cell(11, 1, 2, 7, 12);
  editor.add_cell(12, 4, 7, 8, 12);
  editor.add_cell(13, 4, 7, 9, 12);
  editor.add_cell(14, 5, 6, 10, 12);
  editor.add_cell(15, 5, 6, 11, 12);
  editor.add_cell(16, 8, 11, 0, 12);
  editor.add_cell(17, 8, 11, 1, 12);
  editor.add_cell(18, 9, 10, 2, 12);
  editor.add_cell(19, 9, 10, 3, 12);

  editor.close();

  // Refine base mesh
  for (unsigned int j = 0; j < nrefine; ++j)
  {
    Mesh mesh2(base_mesh.mpi_comm());
    refine(mesh2, base_mesh, false);
    base_mesh = mesh2;
  }

  // Copy over all points from base mesh to mesh

  MeshEditor editor2;
  editor2.open(mesh, tdim, gdim, degree);
  editor2.init_vertices_global(base_mesh.num_vertices(), base_mesh.num_vertices());
  for (VertexIterator v(base_mesh); !v.end(); ++v)
    editor2.add_vertex(v->index(), v->point());
  editor2.init_cells_global(base_mesh.num_cells(), base_mesh.num_cells());
  for (CellIterator c(base_mesh); !c.end(); ++c)
  {
    std::vector<std::size_t> idx(c->entities(0), c->entities(0)+4);
    editor2.add_cell(c->index(), idx);
  }

  // Lift all surface vertices to spherical surface
  move_surface_points(mesh);

  if (degree == 2)
  {
    // Initialise entities required for this degree polynomial mesh
    // and allocate space for the point coordinate data
    editor2.init_entities();

    // Mark all edges which are on the surface
    mesh.init(tdim - 1, tdim);
    std::vector<bool> emarker(mesh.num_edges(), false);
    for (FacetIterator f(mesh); !f.end(); ++f)
    {
      if (f->num_entities(tdim) == 1)
      {
        // Surface facet
        for (EdgeIterator e(*f); !e.end(); ++e)
          emarker[e->index()] = true;
      }
    }

    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      Point v0 = Vertex(mesh, e->entities(0)[0]).point();
      Point pt = e->midpoint();

      // Move any surface nodes to curved surface
      if (emarker[e->index()])
        pt *= v0.norm()/pt.norm();

      // Add Edge-based point
      editor2.add_entity_point(1, 0, e->index(), pt);
    }
  }

  editor2.close();
}
//-----------------------------------------------------------------------------
void UnitSphereMesh::move_surface_points(Mesh& mesh)
{
  const std::size_t tdim = mesh.topology().dim();

  mesh.init(tdim - 1, tdim);

  // Mark all vertices which are on the surface
  std::vector<bool> marker(mesh.num_vertices(), false);
  for (FacetIterator f(mesh); !f.end(); ++f)
  {
    if (f->num_entities(tdim) == 1)
    {
      // Surface facet
      for (VertexIterator v(*f); !v.end(); ++v)
        marker[v->index()] = true;
    }
  }

  const std::size_t gdim = mesh.geometry().dim();
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    std::size_t idx = v->index();
    if (marker[idx])
    {
      Point pt = v->point();
      pt *= 1.0/pt.norm();
      std::copy(pt.coordinates(),
                pt.coordinates() + gdim,
                mesh.geometry().x().begin() + idx*gdim);
    }
  }
}
