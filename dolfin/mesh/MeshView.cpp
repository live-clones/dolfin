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
  mv_index.resize(tdim + 1);
  mv_index[tdim] = indices;

  // Set cell type, tdim, gdim here
  MeshEditor editor;
  editor.open(*this, mvmesh->type().cell_type(), tdim,
              mvmesh->geometry().dim());
  editor.close();

  // Reverse mapping from full Mesh to MeshView
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

  MeshTopology& _topo = topology();
  _topo(tdim, 0).set(new_cells);
  _topo.init(tdim, indices.size(), MPI::sum(mpi_comm(), indices.size()));
  // FIXME: set num_global_vertices correctly (some are shared)
  _topo.init(0, vertex_num, MPI::sum(mpi_comm(), vertex_num));

  // Generate some global indices
  const std::size_t cell_offset
    = MPI::global_offset(mpi_comm(), indices.size(), true);
  _topo.init_global_indices(tdim, indices.size());
  for (std::size_t i = 0; i != indices.size(); ++i)
    _topo.set_global_index(tdim, i, i + cell_offset);

  const std::size_t vertex_offset
    = MPI::global_offset(mpi_comm(), vertex_num, true);
  _topo.init_global_indices(0, vertex_num);
  for (std::size_t i = 0; i != vertex_num; ++i)
    _topo.set_global_index(0, i, i + vertex_offset);

  // FIXME: how to deal with ghosts?
  _topo.init_ghost(0, vertex_num);
  _topo.init_ghost(tdim, indices.size());

  // Copy geometry over
  geometry() = mvmesh->geometry();
  geometry().set_indices(mv_index[0]);
}
//-----------------------------------------------------------------------------
