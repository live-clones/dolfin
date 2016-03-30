// Copyright (C) 2016 Chris Richardson
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

#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/refinement/PlazaRefinementND.h>
#include <dolfin/geometry/Plane.h>

#include "MeshSlicer.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Mesh MeshSlicer::add(const Mesh& mesh, const Plane& surface)
{
  if (mesh.geometry().degree() != 1)
  {
    dolfin_error("MeshSlicer.cpp",
                 "slice higher-order geometry meshes",
                 "Not implemented");
  }

  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t tdim = mesh.topology().dim();

  if (tdim != 3 and gdim != 3)
  {
    dolfin_error("MeshSlicer.cpp",
                 "slice mesh",
                 "Only implemented for slicing 3D meshes");
  }

  // Assign new vertices where surface cuts mesh Edges, and store their
  // coordinates
  std::map<std::size_t, std::size_t> edge_to_vertex;

  std::vector<double> new_geometry(mesh.coordinates().begin(),
                                   mesh.coordinates().end());

  std::size_t count = mesh.num_vertices();
  for (EdgeIterator e(mesh); !e.end(); ++e)
  {
    std::pair<bool, Point> intersect = surface.intersection(*e);
    if (intersect.first)
    {
      const Point& p = intersect.second;
      new_geometry.insert(new_geometry.end(),
                          p.coordinates(), p.coordinates() + gdim);
      edge_to_vertex[e->index()] = count;
      ++count;
    }
  }
  boost::multi_array_ref<double, 2> geometry(new_geometry.data(),
                                             boost::extents
                                             [count][gdim]);

  std::vector<std::size_t> topology;
  std::vector<std::size_t> local_map(10);
  for (CellIterator cell(mesh); !cell.end(); ++cell)
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
        double e_len = Edge(mesh,
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

    std::vector<std::size_t> tetrahedra;
    PlazaRefinementND::get_simplices(tetrahedra, marked_edges,
                                     facet_key_edge, 3, false);

    // Add to topology, mapping to the new local indexing
    for (auto p : tetrahedra)
      topology.push_back(local_map[p]);
  }

  Mesh _mesh;
  create_mesh(_mesh, geometry, topology, tdim);
  return _mesh;
}
//-----------------------------------------------------------------------------
Mesh MeshSlicer::cut(const Mesh& mesh, const Plane& surface)
{

  if (mesh.geometry().degree() != 1)
  {
    dolfin_error("MeshSlicer.cpp",
                 "slice higher-order geometry meshes",
                 "Not implemented");
  }

  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t tdim = mesh.topology().dim();

  if (tdim != 3 and gdim != 3)
  {
    dolfin_error("MeshSlicer.cpp",
                 "slice mesh",
                 "Only implemented for slicing 3D meshes");
  }

  // Assign new vertices where surface cuts mesh Edges, and store their
  // coordinates
  std::map<std::size_t, std::size_t> edge_to_vertex;
  std::vector<double> new_geometry;
  std::size_t count = 0;
  for (EdgeIterator e(mesh); !e.end(); ++e)
  {
    std::pair<bool, Point> intersect = surface.intersection(*e);
    if (intersect.first)
    {
      const Point& p = intersect.second;
      new_geometry.insert(new_geometry.end(),
                          p.coordinates(), p.coordinates() + gdim);
      edge_to_vertex[e->index()] = count;
      ++count;
    }
  }

  mesh.init(tdim, 1);
  std::vector<std::size_t> topology;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    std::vector<size_t> local_topology;
    for (unsigned int i = 0; i != 6; ++i)
    {
      const std::size_t edge_idx = cell->entities(1)[i];
      auto map_it = edge_to_vertex.find(edge_idx);
      if (map_it != edge_to_vertex.end())
        local_topology.push_back(map_it->second);
    }

    if (local_topology.size() == 3)
    {
      topology.insert(topology.end(),
                      local_topology.begin(),
                      local_topology.end());
    }
    else if(local_topology.size() == 4)
    {
      topology.insert(topology.end(),
                      local_topology.begin(),
                      local_topology.begin() + 3);
      topology.insert(topology.end(),
                      local_topology.begin() + 1,
                      local_topology.begin() + 4);
    }

  }

  boost::multi_array_ref<double, 2> geometry(new_geometry.data(),
                                             boost::extents
                                             [count][gdim]);

  Mesh _mesh;
  create_mesh(_mesh, geometry, topology, tdim - 1);
  return _mesh;
}
//-----------------------------------------------------------------------------
void MeshSlicer::create_mesh(Mesh& mesh,
                             const boost::multi_array<double, 2>& geometry,
                             const std::vector<std::size_t>& topology,
                             std::size_t tdim)
{
  MeshEditor ed;

  // Create a mesh of triangles or tetrahedra in 3D
  ed.open(mesh, tdim, 3);

  ed.init_vertices(geometry.size());
  unsigned int i = 0;
  for (const auto p : geometry)
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
