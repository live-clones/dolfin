// Copyright (C) 2015 Chris Richardson
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

#include "MeshView.h"
#include "MeshEditor.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshView::MeshView(std::shared_ptr<Mesh> mvmesh, std::size_t tdim,
         const std::vector<std::size_t>& indices)
  : Mesh(mvmesh->mpi_comm()), mv_mesh(mvmesh)
{
  dolfin_assert(tdim <= mvmesh->topology().dim());

  mv_index.resize(tdim + 1);
  mv_index[tdim] = indices;

  // Set cell type, tdim, gdim here
  MeshEditor editor;
  editor.open(*this, mvmesh->type().entity_type(tdim), tdim,
              mvmesh->geometry().dim());
  editor.close();

  // Reverse mapping from full Mesh to MeshView for vertices
  std::map<std::size_t, std::size_t> mv_index_map;
  const MeshConnectivity& mvconn = mvmesh->topology()(tdim, 0);

  std::size_t vertex_num = 0;
  std::vector<std::vector<std::size_t> > new_cells(indices.size());
  for (unsigned int j = 0; j != indices.size(); ++j)
  {
    const std::size_t idx = indices[j];
    const std::size_t num_cell_vertices = mvconn.size(idx);
    for (unsigned int i = 0; i != num_cell_vertices; ++i)
    {
      const std::size_t main_idx = mvconn(idx)[i];
      auto mapit = mv_index_map.insert(std::pair<std::size_t, std::size_t>
                                       (main_idx, vertex_num));
      if (mapit.second)
      {
        // Insert forward mapping for vertices
        mv_index[0].push_back(main_idx);
        ++vertex_num;
      }
      new_cells[j].push_back(mapit.first->second);
    }
  }

  // Sort out shared vertices
  const std::size_t mpi_size = MPI::size(mpi_comm());
  const auto& mvshared = mvmesh->topology().shared_entities(0);
  auto& view_shared = topology().shared_entities(0);
  const auto& mvglobal = mvmesh->topology().global_indices(0);
  const std::size_t num_global_main = mvmesh->topology().size_global(0);
  std::size_t local_count = 0;
  std::vector<std::size_t> vertex_global_index(vertex_num);

  // For MPI::all_to_all to communicate shared vertex numbering
  std::vector<std::vector<std::size_t>> send_vertex_numbering(mpi_size);
  std::vector<std::vector<std::size_t>> recv_vertex_numbering(mpi_size);

  // Mapping back from main global index to local view index
  std::map<std::size_t, std::size_t> main_global_to_view;

  for (unsigned int i = 0; i != vertex_num; ++i)
  {
    const std::size_t main_idx = mv_index[0][i];
    const auto shared_it = mvshared.find(main_idx);
    if (shared_it != mvshared.end())
    {
      view_shared.insert
        (std::pair<unsigned int, std::set<unsigned int>>(i, shared_it->second));
      // Send to index owner, for consistent numbering across processes
      const std::size_t dest = MPI::index_owner(mpi_comm(), mvglobal[main_idx], num_global_main);
      send_vertex_numbering[dest].push_back(mvglobal[main_idx]);
      main_global_to_view.insert(std::pair<std::size_t, std::size_t>
                                 (mvglobal[main_idx], i));
    }
    else
    {
      // Unshared vertex, numbered locally
      vertex_global_index[i] = local_count;
      ++local_count;
    }
  }

  // Communicate remotely numbered indices
  MPI::all_to_all(mpi_comm(), send_vertex_numbering, recv_vertex_numbering);
  std::map<std::size_t, std::size_t> global_index_set;
  // Create map from all received main global indices to local_count
  for (auto &proc : recv_vertex_numbering)
    for (auto &q : proc)
    {
      auto map_insert = global_index_set.insert(std::pair<std::size_t, std::size_t>(q, local_count));
      if (map_insert.second)
        ++local_count;
    }

  // Correct for global vertex index offset
  const std::size_t vertex_offset
    = MPI::global_offset(mpi_comm(), local_count, true);
  for (auto &v : vertex_global_index)
    v += vertex_offset;
  for (auto &v : global_index_set)
    v.second += vertex_offset;

  // Convert incoming main global index to meshview global index
  for (auto &proc : recv_vertex_numbering)
    for (auto &q : proc)
    {
      auto map_it = global_index_set.find(q);
      dolfin_assert(map_it != global_index_set.end());
      q = map_it->second;
    }

  // Send reply back to originator
  std::vector<std::vector<std::size_t>> reply_vertex_numbering(mpi_size);
  MPI::all_to_all(mpi_comm(), recv_vertex_numbering, reply_vertex_numbering);

  // Convert main global back to meshview local and save new meshview global
  // index
  for (unsigned int i = 0; i != mpi_size; ++i)
  {
    std::vector<std::size_t>& send_v = send_vertex_numbering[i];
    std::vector<std::size_t>& reply_v = reply_vertex_numbering[i];
    dolfin_assert(send_v.size() == reply_v.size());

    for (unsigned int j = 0; j != send_v.size(); ++j)
    {
      auto map_it = main_global_to_view.find(send_v[j]);
      dolfin_assert(map_it != main_global_to_view.end());
      vertex_global_index[map_it->second]
        = reply_v[j];
    }
  }

  // FIXME: some cells are shared if ghosts are used
  MeshTopology& _topo = topology();
  _topo(tdim, 0).set(new_cells);
  _topo.init(tdim, indices.size(), MPI::sum(mpi_comm(), indices.size()));

  // Generate some global indices for cells
  const std::size_t cell_offset
    = MPI::global_offset(mpi_comm(), indices.size(), true);
  _topo.init_global_indices(tdim, indices.size());
  for (std::size_t i = 0; i != indices.size(); ++i)
    _topo.set_global_index(tdim, i, i + cell_offset);

  const std::size_t num_global_vertices = MPI::sum(mpi_comm(), local_count);
  _topo.init(0, vertex_num, num_global_vertices);
  _topo.init_global_indices(0, vertex_num);
  for (std::size_t i = 0; i != vertex_num; ++i)
    _topo.set_global_index(0, i, vertex_global_index[i]);

  // FIXME: how to deal with ghosts?
  _topo.init_ghost(0, vertex_num);
  _topo.init_ghost(tdim, indices.size());

  // Copy geometry over
  geometry() = mvmesh->geometry();
  geometry().set_indices(mv_index[0]);
}
//-----------------------------------------------------------------------------
