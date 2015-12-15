// Copyright (C) 2015 Chris Richardson
//
// This file is part of Narwal
//
// Narwal is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Narwal is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Narwal. If not, see <http://www.gnu.org/licenses/>.

#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/refinement/PlazaRefinementND.h>
#include <dolfin/geometry/Plane.h>

#include "SubTriangulation.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
SubTriangulation::SubTriangulation(const Mesh& mesh,
                                   const Plane& surface)
  : _mesh(mesh)
{
  dolfin_assert(mesh.cell_type() == CellType::tetrahedron);
  build(surface);
}
//-----------------------------------------------------------------------------
SubTriangulation::SubTriangulation(const Cell& cell,
                                   const Plane& surface)
{
  dolfin_assert(cell.type() == CellType::tetrahedron);

  boost::multi_array<double, 2> geometry (boost::extents[4][3]);
  std::vector<std::size_t> topology = {0, 1, 2, 3};

  unsigned int i = 0;
  for (VertexIterator v(cell); !v.end(); ++v)
  {
    std::copy(v->x(), v->x() + 3, geometry[i].begin());
    ++i;
  }
  create_mesh(_mesh, geometry, topology, 3);

  build(surface);
}
//-----------------------------------------------------------------------------
void SubTriangulation::build(const Plane& surface)
{
  // Assign new vertices where surface cuts mesh Edges, and store their
  // coordinates
  std::map<std::size_t, std::size_t> edge_to_vertex;

  std::vector<double> new_geometry(_mesh.coordinates().begin(),
                                   _mesh.coordinates().end());
  std::size_t count = _mesh.num_vertices();
  _mesh.init(1, 0);

  for (EdgeIterator e(_mesh); !e.end(); ++e)
  {
    std::pair<bool, Point> intersect = surface.intersection(*e);
    if (intersect.first)
    {
      const Point& p = intersect.second;
      new_geometry.insert(new_geometry.end(),
                          p.coordinates(), p.coordinates() + 3);
      edge_to_vertex[e->index()] = count;
      ++count;
    }
  }

  boost::multi_array_ref<double, 2> geometry(new_geometry.data(),
                                             boost::extents
                                             [count][3]);

  std::vector<std::size_t> topology;
  std::vector<std::size_t> local_map(10);
  for (CellIterator cell(_mesh); !cell.end(); ++cell)
  {
    // Create local mapping for this cell - add existing vertices
    for (unsigned int i = 0; i != 4; ++i)
      local_map[i] = cell->entities(0)[i];

    // Also add vertices where edges are cut by surface
    std::vector<bool> marked_edges(6);
    for (unsigned int i = 0; i != 6; ++i)
    {
      const std::size_t edge_idx = cell->entities(1)[i];
      auto map_it = edge_to_vertex.find(edge_idx);
      if(map_it != edge_to_vertex.end())
      {
        local_map[4 + i] = map_it->second;
        marked_edges[i] = true;
      }
    }

    // Edges on each facet in cell local indexing
    static std::size_t facet_edge[4][3]
      = { {1, 2, 0},
          {3, 4, 0},
          {3, 5, 1},
          {4, 5, 2} };

    std::vector<std::size_t> facet_key_edge;
    for (std::size_t i = 0; i != 4; ++i)
    {
      std::size_t max_idx = 0;
      double max_len = 0.0;
      bool max_mark = false;

      for (unsigned int j = 0; j != 3; ++j)
      {
        std::size_t e_idx = facet_edge[i][j];
        double e_len = Edge(_mesh,
                            cell->entities(1)[e_idx]).length();
        bool e_mark = marked_edges[e_idx];

        if (e_mark > max_mark)
        {
          max_mark = e_mark;
          max_len = e_len;
          max_idx = e_idx;
        }
        else if (e_mark == max_mark && e_len > max_len)
        {
          max_len = e_len;
          max_idx = e_idx;
        }
      }
      facet_key_edge.push_back(max_idx);
    }

    std::vector<std::vector<std::size_t> > tetrahedra;
    PlazaRefinementND::get_simplices(tetrahedra, marked_edges,
                                     facet_key_edge, 3);

    // Add to topology, mapping to the new local indexing
    for (auto tet : tetrahedra)
      for (auto p : tet)
        topology.push_back(local_map[p]);

  }

  create_mesh(_mesh, geometry, topology, 3);
}
//-----------------------------------------------------------------------------
void SubTriangulation::create_mesh(
  Mesh& mesh,
  const boost::multi_array<double, 2>& geometry,
  const std::vector<std::size_t>& topology,
  std::size_t tdim)
{
  MeshEditor ed;

  // Create a mesh of triangles
  ed.open(mesh, tdim, 3);

  ed.init_vertices(geometry.size());
  unsigned int i = 0;
  for (auto p : geometry)
  {
    ed.add_vertex(i, p[0], p[1], p[2]);
    ++i;
  }

  dolfin_assert(topology.size()%(tdim + 1) == 0);
  ed.init_cells(topology.size()/(tdim + 1));
  i = 0;
  for (auto t = topology.begin(); t != topology.end(); t += (tdim + 1))
  {
    if (tdim == 3)
      ed.add_cell(i, *t, *(t+1), *(t+2), *(t+3));
    else if (tdim == 2)
      ed.add_cell(i, *t, *(t+1), *(t+2));
    ++i;
  }

  ed.close();
}
//-----------------------------------------------------------------------------
