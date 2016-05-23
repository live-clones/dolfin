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
//

#include "Facet.h"
#include "Mesh.h"
#include "MeshEditor.h"
#include "MeshRelation.h"

#include "NewBoundaryMesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::shared_ptr<const MeshRelation>
NewBoundaryMesh::create_boundary(std::shared_ptr<const Mesh> mesh)
{
  const unsigned int tdim = mesh->topology().dim();

  // Make sure facet-cell connection is initialised
  mesh->init(tdim - 1, tdim);

  // Get indices of surface facets
  std::vector<std::size_t> indices;
  for (FacetIterator f(*mesh); !f.end(); ++f)
  {
    if (f->num_global_entities(tdim) == 1)
      indices.push_back(f->index());
  }

  return create(mesh, indices, tdim - 1);
}
//-----------------------------------------------------------------------------
std::shared_ptr<const MeshRelation>
NewBoundaryMesh::create(std::shared_ptr<const Mesh> mesh,
                        std::vector<std::size_t> indices,
                        std::size_t tdim)
{
  std::cout << "tdim = " << tdim << "\n";
  std::cout << "num indices = " << indices.size() << "\n";

  auto boundary = std::make_shared<Mesh>();

  MeshEditor editor;
  editor.open(*boundary, mesh->type().entity_type(tdim), tdim,
              mesh->geometry().dim());

  // Get number of cell vertices on boundary cell
  std::size_t num_cell_vertices = mesh->type().num_vertices(tdim);

  // Reverse mapping from full Mesh for vertices
  std::map<std::size_t, std::size_t> vertex_rev_map;
  // Forward map to full Mesh for vertices
  std::vector<std::size_t> vertex_fwd_map;

  // Establish a new local numbering for all vertices
  // in "boundary" and the mapping back to the local
  // numbering in the main mesh

  std::size_t vertex_num = 0;
  // FIXME: may fail with ghost mesh - assumes no shared cells
  editor.init_cells_global(indices.size(),
                           MPI::sum(mesh->mpi_comm(), indices.size()));
  const auto& conn_tdim = mesh->topology()(tdim, 0);
  std::vector<std::size_t> new_cell;
  for (unsigned int j = 0; j != indices.size(); ++j)
  {
    const std::size_t idx = indices[j];
    new_cell.clear();
    for (unsigned int i = 0; i != num_cell_vertices; ++i)
    {
      const std::size_t main_idx = conn_tdim(idx)[i];
      auto mapit = vertex_rev_map.insert({main_idx, vertex_num});
      if (mapit.second)
      {
        vertex_fwd_map.push_back(main_idx);
        ++vertex_num;
      }
      new_cell.push_back(mapit.first->second);
    }
    std::cout << "j = " << j << " ";
    for (auto q : new_cell)
      std::cout << q <<" ";
    std::cout <<" \n";

    editor.add_cell(j, new_cell);
  }

  // Initialise vertex global indices on boundary

  const std::size_t mpi_size = MPI::size(mesh->mpi_comm());
  const std::size_t mpi_rank = MPI::rank(mesh->mpi_comm());

  // Create sharing map for boundary from main mesh data
  auto& bound_shared = boundary->topology().shared_entities(0);
  const auto& mesh_shared = mesh->topology().shared_entities(0);
  const auto& mesh_global = mesh->topology().global_indices(0);

  // Global numbering on boundary
  std::size_t local_count = 0;
  std::vector<std::size_t> vertex_global_index(vertex_num);
  std::vector<std::vector<std::size_t>> send_vertex_numbering(mpi_size);
  std::vector<std::vector<std::size_t>> recv_vertex_numbering(mpi_size);

  // Map from global vertices on main mesh to local on Boundary
  std::map<std::size_t, std::size_t> main_global_to_bound;

  for (unsigned int i = 0; i != vertex_num; ++i)
  {
    const std::size_t main_idx = vertex_fwd_map[i];
    const auto shared_it = mesh_shared.find(main_idx);
    if (shared_it == mesh_shared.end())
    {
      // Vertex not shared in main mesh, so number locally
      vertex_global_index[i] = local_count;
      ++local_count;
    }
    else
    {
      bound_shared.insert({i, shared_it->second});
      // Send to remote 'owner' for numbering.
      std::size_t dest = *(shared_it->second.begin());
      if (dest > mpi_rank)
      {
        // Shared, but local - number locally
        vertex_global_index[i] = local_count;
        ++local_count;
      }
      else
        send_vertex_numbering[dest].push_back(mesh_global[main_idx]);

      main_global_to_bound.insert({mesh_global[main_idx], i});
    }
  }

  // Send global indices of all unnumbered vertices to owner
  MPI::all_to_all(mesh->mpi_comm(), send_vertex_numbering, recv_vertex_numbering);

  // Create global->local map for *all* shared vertices of main mesh
  std::map<std::size_t, std::size_t> main_global_to_local;
  for (auto it : mesh_shared)
    main_global_to_local.insert({mesh_global[it.first], it.first});

  // Search for received vertices which are not already there. This may be because
  // they are not part of a Boundary cell for this process.
  for (auto &p : recv_vertex_numbering)
  {
    for (auto &q : p)
    {
      if (main_global_to_bound.find(q) == main_global_to_bound.end())
      {
        // Not found - shared, but not part of local MeshView, yet.
        // This Vertex may have no locally associated Cell on the boundary
        auto mgl_find = main_global_to_local.find(q);
        dolfin_assert(mgl_find != main_global_to_local.end());
        const std::size_t main_idx = mgl_find->second;
        vertex_fwd_map.push_back(main_idx);
        vertex_global_index.push_back(local_count);
        main_global_to_bound.insert({q, vertex_num});
        const auto shared_it = mesh_shared.find(main_idx);
        dolfin_assert(shared_it != mesh_shared.end());
        bound_shared.insert({vertex_num, shared_it->second});
        ++local_count;
        ++vertex_num;
      }
    }
  }

  // All shared vertices are now numbered

  // Add geometry and close "boundary" mesh
  editor.init_vertices_global(vertex_fwd_map.size(),
                              MPI::sum(mesh->mpi_comm(), local_count));
  for (unsigned int i = 0; i != vertex_fwd_map.size(); ++i)
    editor.add_vertex(i, mesh->geometry().point(vertex_fwd_map[i]));
  editor.close();

  // Correct for global vertex index offset
  const std::size_t vertex_offset
    = MPI::global_offset(mesh->mpi_comm(), local_count, true);
  for (auto &v : vertex_global_index)
    v += vertex_offset;

  // Convert incoming main global index to Boundary global index
  for (auto &p : recv_vertex_numbering)
    for (auto &q : p)
    {
      const auto map_it = main_global_to_bound.find(q);
      dolfin_assert(map_it != main_global_to_bound.end());
      q = vertex_global_index[map_it->second];
    }

  // Send reply back to originator
  std::vector<std::vector<std::size_t>> reply_vertex_numbering(mpi_size);
  MPI::all_to_all(mesh->mpi_comm(), recv_vertex_numbering, reply_vertex_numbering);

  // Convert main global back to Boundary local and save new Boundary global
  // index
  for (unsigned int i = 0; i != mpi_size; ++i)
  {
    std::vector<std::size_t>& send_v = send_vertex_numbering[i];
    std::vector<std::size_t>& reply_v = reply_vertex_numbering[i];
    dolfin_assert(send_v.size() == reply_v.size());

    for (unsigned int j = 0; j != send_v.size(); ++j)
    {
      auto map_it = main_global_to_bound.find(send_v[j]);
      dolfin_assert(map_it != main_global_to_bound.end());
      vertex_global_index[map_it->second] = reply_v[j];
    }
  }

  MeshTopology& btopo = boundary->topology();

  // Generate some global indices for cells
  // FIXME: some cells may be shared if ghosts are used
  const std::size_t cell_offset
    = MPI::global_offset(mesh->mpi_comm(), indices.size(), true);
  btopo.init_global_indices(tdim, indices.size());
  for (std::size_t i = 0; i != indices.size(); ++i)
    btopo.set_global_index(tdim, i, i + cell_offset);

  // Set global vertex indices
  btopo.init_global_indices(0, vertex_num);
  for (std::size_t i = 0; i != vertex_num; ++i)
    btopo.set_global_index(0, i, vertex_global_index[i]);

  // Record relation between meshes
  auto relation = std::make_shared<MeshRelation>(mesh, boundary);
  //  relation->init(0, vertex_fwd_map);
  //  relation->init(tdim, indices);

  return relation;
}
//-----------------------------------------------------------------------------
