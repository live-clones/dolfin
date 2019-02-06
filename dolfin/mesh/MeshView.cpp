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
#include "Mesh.h"
#include "MeshEditor.h"
#include "MeshFunction.h"
#include "MeshEntityIterator.h"

#include "MeshView.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Mesh MeshView::create(const MeshFunction<std::size_t>& marker,
                      std::size_t tag)
{
  // Get original Mesh and tdim of marker
  std::shared_ptr<const Mesh> mesh = marker.mesh();
  unsigned int tdim = marker.dim();

  // Get indices of marked entities - each of these will represent a Cell in "new_mesh"
  std::vector<std::size_t> indices;
  for (std::size_t idx = 0; idx < mesh->num_entities(tdim); ++idx)
  {
    if (marker[idx] == tag)
      indices.push_back(idx);
  }

  // Create a new Mesh of dimension tdim
  Mesh new_mesh;
  MeshEditor editor;
  editor.open(new_mesh, mesh->type().entity_type(tdim), tdim,
              mesh->geometry().dim());

  std::size_t num_cell_vertices = mesh->type().num_vertices(tdim);

  // Reverse mapping from full Mesh for vertices
  std::map<std::size_t, std::size_t> vertex_rev_map;
  // Forward map to full Mesh for vertices
  std::vector<std::size_t> vertex_fwd_map;

  // Establish a new local numbering for all vertices in "new_mesh"
  // and create the mapping back to the local numbering in "mesh"

  std::size_t vertex_num = 0;
  // FIXME: may fail with ghost mesh - assumes no shared cells
  editor.init_cells_global(indices.size(),
                           MPI::sum(mesh->mpi_comm(), indices.size()));
  const auto& conn_tdim = mesh->topology()(tdim, 0);

  // Add cells to new_mesh
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
    editor.add_cell(j, new_cell);
  }

  // Initialise vertex global indices in new_mesh
  const std::size_t mpi_size = MPI::size(mesh->mpi_comm());
  const std::size_t mpi_rank = MPI::rank(mesh->mpi_comm());

  // Create sharing map from main mesh data
  auto& new_shared = new_mesh.topology().shared_entities(0);
  const auto& mesh_shared = mesh->topology().shared_entities(0);
  const auto& mesh_global = mesh->topology().global_indices(0);

  // Global numbering in new_mesh. Shared vertices are numbered by
  // the lowest rank sharing process
  std::size_t local_count = 0;
  std::vector<std::size_t> vertex_global_index(vertex_num);
  std::vector<std::vector<std::size_t>> send_vertex_numbering(mpi_size);
  std::vector<std::vector<std::size_t>> recv_vertex_numbering(mpi_size);

  // Map from global vertices on main mesh to local on new_mesh
  std::map<std::size_t, std::size_t> main_global_to_new;

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
      // Dof is owned by the lowest rank process
      //std::size_t dest = *std::min_element(shared_it->second.begin(), shared_it->second.end());
      std::size_t dest = *(shared_it->second.begin());
      new_shared.insert({i, shared_it->second});

      if (dest > mpi_rank)
      {
        // Shared, but local - number locally
        vertex_global_index[i] = local_count;
        ++local_count;
      }
      else
        send_vertex_numbering[dest].push_back(mesh_global[main_idx]);

      main_global_to_new.insert({mesh_global[main_idx], i});
    }
  }

  // Send global indices of all unnumbered vertices to owner
  MPI::all_to_all(mesh->mpi_comm(), send_vertex_numbering, recv_vertex_numbering);

  // Create global->local map for *all* shared vertices of main mesh
  std::map<std::size_t, std::size_t> main_global_to_local;
  for (auto it : mesh_shared)
    main_global_to_local.insert({mesh_global[it.first], it.first});

  // Search for received vertices which are not already there. This may be because
  // they are not part of a new_mesh cell for this process.
  for (auto &p : recv_vertex_numbering)
  {
    for (auto &q : p)
    {
      if (main_global_to_new.find(q) == main_global_to_new.end())
      {
        // Not found - shared, but not part of local MeshView, yet.
        // This Vertex may have no locally associated Cell in new_mesh
        auto mgl_find = main_global_to_local.find(q);
        dolfin_assert(mgl_find != main_global_to_local.end());
        const std::size_t main_idx = mgl_find->second;
        vertex_fwd_map.push_back(main_idx);
        vertex_global_index.push_back(local_count);
        main_global_to_new.insert({q, vertex_num});
        const auto shared_it = mesh_shared.find(main_idx);
        dolfin_assert(shared_it != mesh_shared.end());
        new_shared.insert({vertex_num, shared_it->second});
#if 0
        // DEBUG : Print new shared vertex global index
        for(auto ii=shared_it->second.begin(); ii!=shared_it->second.end(); ++ii)
          std::cout << "new_shared - "
                    << " new mesh index = " << vertex_num
                    << ", vertex_global_index = " << local_count
                    << ", main global index = " << q
                    << ", owned by " << *ii << std::endl;
#endif
        ++local_count;
        ++vertex_num;
      }
    }
  }

  // All shared vertices are now numbered

  // Add geometry and close new_mesh
  editor.init_vertices_global(vertex_fwd_map.size(),
                              MPI::sum(mesh->mpi_comm(), local_count));
  for (unsigned int i = 0; i != vertex_fwd_map.size(); ++i)
    new_mesh.geometry().set(i, mesh->geometry().point(vertex_fwd_map[i]).coordinates());

  editor.close();

  // Correct for global vertex index offset
  const std::size_t vertex_offset
    = MPI::global_offset(mesh->mpi_comm(), local_count, true);
  for (auto &v : vertex_global_index)
    v += vertex_offset;

  // Convert incoming main global index to new_mesh global index
  for (auto &p : recv_vertex_numbering)
    for (auto &q : p)
    {
      const auto map_it = main_global_to_new.find(q);
      dolfin_assert(map_it != main_global_to_new.end());
      q = vertex_global_index[map_it->second];
    }

  // Send reply back to originator
  std::vector<std::vector<std::size_t>> reply_vertex_numbering(mpi_size);
  MPI::all_to_all(mesh->mpi_comm(), recv_vertex_numbering, reply_vertex_numbering);

  // Convert main global back to new_mesh local and save new_mesh global
  // index
  for (unsigned int i = 0; i != mpi_size; ++i)
  {
    std::vector<std::size_t>& send_v = send_vertex_numbering[i];
    std::vector<std::size_t>& reply_v = reply_vertex_numbering[i];
    dolfin_assert(send_v.size() == reply_v.size());

    for (unsigned int j = 0; j != send_v.size(); ++j)
    {
      auto map_it = main_global_to_new.find(send_v[j]);
      dolfin_assert(map_it != main_global_to_new.end());
      vertex_global_index[map_it->second] = reply_v[j];
    }
  }

  MeshTopology& new_topo = new_mesh.topology();

  // Generate some global indices for cells
  // FIXME: some cells may be shared if ghosts are used
  const std::size_t cell_offset
    = MPI::global_offset(mesh->mpi_comm(), indices.size(), true);
  new_topo.init_global_indices(tdim, indices.size());
  for (std::size_t i = 0; i != indices.size(); ++i)
    new_topo.set_global_index(tdim, i, i + cell_offset);

  // Set global vertex indices
  new_topo.init_global_indices(0, vertex_num);
  for (std::size_t i = 0; i != vertex_num; ++i)
    new_topo.set_global_index(0, i, vertex_global_index[i]);

  // Store relationship between meshes
  new_topo._mapping.insert(std::make_pair(mesh->id(), std::make_shared<MeshView>(mesh, vertex_fwd_map, indices)));

  return new_mesh;
}
