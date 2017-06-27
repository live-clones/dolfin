// Copyright (C) 2017 Nate Sime and Chris Richardson
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

#include <dolfin/common/ArrayView.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>

#include "BoundingBoxTree.h"
#include "CollisionDetection.h"
#include "Point.h"

#include "GeometricContact.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::vector<Point> GeometricContact::create_deformed_segment_volume_3d(Mesh& mesh, std::size_t facet_index, const Function& u)
{
  Facet facet(mesh, facet_index);
  Point X1 = Vertex(mesh, facet.entities(0)[0]).point();
  Point X2 = Vertex(mesh, facet.entities(0)[1]).point();
  Point X3 = Vertex(mesh, facet.entities(0)[2]).point();

  // Get id of attached cell
  std::size_t id = facet.entities(mesh.topology().dim())[0];

  Array<double> uval(3);

  const Cell cell(mesh, id);
  ufc::cell ufc_cell;
  cell.get_cell_data(ufc_cell);

  const Array<double> _x1(3, X1.coordinates());
  u.eval(uval, _x1, cell, ufc_cell);
  Point x1 = X1 + Point(uval);

  const Array<double> _x2(3, X2.coordinates());
  u.eval(uval, _x2, cell, ufc_cell);
  Point x2 = X2 + Point(uval);

  const Array<double> _x3(3, X3.coordinates());
  u.eval(uval, _x3, cell, ufc_cell);
  Point x3 = X3 + Point(uval);

  return {X1, X2, X3, x1, x2, x3};
}
//-----------------------------------------------------------------------------
bool GeometricContact::check_tet_set_collision(const Mesh& master_mesh, std::size_t mindex,
                                               const Mesh& slave_mesh, std::size_t sindex)
{

  for (unsigned int i = mindex; i < mindex + 3; ++i)
    for (unsigned int j = sindex; j < sindex + 3; ++j)
    {
      // FIXME: needs to work with zero-volume tetrahedra
      if (CollisionDetection::collides_tetrahedron_tetrahedron(Cell(master_mesh, i),
                                                               Cell(slave_mesh, j)))
        return true;
    }

  return false;
}
//-----------------------------------------------------------------------------
void GeometricContact::contact_surface_map_volume_sweep_3d(Mesh& mesh, Function& u,
const std::vector<std::size_t>& master_facets, const std::vector<std::size_t>& slave_facets)
{
  // Construct a dictionary mapping master facets to their collided slave counterparts.

  // This algorithm find which slave facets are contained in the volume swept out by
  // the master surface and its displacement.

  // :param mesh: The mesh
  // :param u: The (3D) solution vector field
  // :param master_facets: list of master facet indices
  // :param slave_facets: list of slave facet indices
  // :return: mapping master facet index to slave facet indices

  if (mesh.topology().dim() != 3)
  {
    dolfin_error("GeometricContact.cpp",
                 "find contact surface",
                 "Only implemented in 3D");
  }

  // Ensure BBT built
  auto mesh_bb = mesh.bounding_box_tree();

  // Make sure facet->cell connections are made
  mesh.init(mesh.topology().dim() - 1, mesh.topology().dim());

  // Local mesh of master 'prisms', three tetrahedra created per facet
  Mesh master_mesh(mesh.mpi_comm());
  MeshEditor master_ed;
  master_ed.open(master_mesh, mesh.topology().dim(), mesh.geometry().dim());
  std::size_t nc_local = master_facets.size()*3;
  std::size_t nc_global = MPI::sum(mesh.mpi_comm(), nc_local);
  master_ed.init_cells_global(nc_local, nc_global);
  master_ed.init_vertices_global(nc_local*2, nc_global*2);
  int c = 0;
  int v = 0;

  for (auto &mf : master_facets)
  {
    // Get six points defining prism
    std::vector<Point> master_tet_set = create_deformed_segment_volume_3d(mesh, mf, u);
    // Add three tetrahedra
    for (int i = 0; i < 3; ++i)
      master_ed.add_cell(c+i, v+i, v+i+1, v+i+2, v+i+3);
    c += 3;
    // Add six points
    for (int i = 0; i < 6; ++i)
      master_ed.add_vertex(v+i, master_tet_set[i]);
    v += 6;
  }
  master_ed.close();
  auto master_bb = master_mesh.bounding_box_tree();

  // Local mesh of slave 'prisms', three tetrahedra created per facet
  Mesh slave_mesh(mesh.mpi_comm());
  MeshEditor slave_ed;
  slave_ed.open(slave_mesh, mesh.topology().dim(), mesh.geometry().dim());
  nc_local = slave_facets.size()*3;
  nc_global = MPI::sum(mesh.mpi_comm(), nc_local);
  slave_ed.init_cells_global(nc_local, nc_global);
  slave_ed.init_vertices_global(nc_local*2, nc_global*2);
  c = 0;
  v = 0;
  for (auto &sf : slave_facets)
  {
    std::vector<Point> slave_tet_set = create_deformed_segment_volume_3d(mesh, sf, u);
    for (int i = 0; i < 3; ++i)
      slave_ed.add_cell(c+i, v+i, v+i+1, v+i+2, v+i+3);
    c += 3;
    for (int i = 0; i < 6; ++i)
      slave_ed.add_vertex(v+i, slave_tet_set[i]);
    v += 6;
  }
  slave_ed.close();
  auto slave_bb = slave_mesh.bounding_box_tree();

  _master_to_slave.clear();

  std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());
  std::size_t mpi_size = MPI::size(mesh.mpi_comm());

  // Check each master 'prism' against each slave 'prism'
  // Map is stored as local_master_facet -> [mpi_rank, local_index, mpi_rank, local_index ...]
  // First check locally
  for (unsigned int i = 0; i < master_facets.size(); ++i)
    for (unsigned int j = 0; j < slave_facets.size(); ++j)
    {
      if (check_tet_set_collision(master_mesh, i*3, slave_mesh, j*3))
      {
        const std::size_t mf = master_facets[i];
        const std::size_t sf = slave_facets[j];
        _master_to_slave[mf].push_back(mpi_rank);
        _master_to_slave[mf].push_back(sf);
      }
    }

  // Find which [master global/slave entity] BBs overlap in parallel
  if (mpi_size > 1)
  {
    // Find which master processes collide with which local slave cells
    auto overlap = master_bb->compute_process_entity_collisions(*slave_bb);
    auto master_procs = overlap.first;
    auto slave_cells = overlap.second;

    // Get slave facet indices to send
    // The slave mesh consists of repeated units of three tetrahedra,
    // making the prisms which are projected forward. The "prism index"
    // can be obtained by integer division of the cell index.
    // Each prism corresponds to a slave facet of the original mesh
    std::vector<std::vector<std::size_t>> send_facets(mpi_size);
    for (unsigned int i = 0; i < master_procs.size(); ++i)
    {
      const unsigned int master_rank = master_procs[i];
      // Ignore local (already done)
      if (master_rank != mpi_rank)
      {
        // Get facet from cell index (three tets per prism)
        unsigned int facet = slave_cells[i]/3;
        send_facets[master_rank].push_back(facet);
      }
    }

    // Get unique set of facets to send to each process
    for (unsigned int p = 0; p != mpi_size; ++p)
    {
      std::vector<std::size_t>& v = send_facets[p];
      std::sort(v.begin(), v.end());
      auto last = std::unique(v.begin(), v.end());
      v.erase(last, v.end());
    }

    // Get coordinates of prism (18 doubles in 3D)
    // and convert index of slave mesh back to index of main mesh
    std::vector<std::vector<double>> send_coordinates(mpi_size);
    for (unsigned int p = 0; p != mpi_size; ++p)
    {
      std::vector<double>& coords_p = send_coordinates[p];
      for (auto& q : send_facets[p])
      {
        const double *coord_vals = slave_mesh.geometry().vertex_coordinates(q*6);
        coords_p.insert(coords_p.end(), coord_vals, coord_vals + 18);
        // Convert back to main mesh indexing
        q = slave_facets[q];
      }
    }

    std::vector<std::vector<std::size_t>> recv_facets(mpi_size);
    MPI::all_to_all(mesh.mpi_comm(), send_facets, recv_facets);

    std::vector<std::vector<double>> recv_coordinates(mpi_size);
    MPI::all_to_all(mesh.mpi_comm(), send_coordinates, recv_coordinates);

    for (unsigned int proc = 0; proc != mpi_size; ++proc)
    {
      auto& rfacet = recv_facets[proc];
      auto& coord = recv_coordinates[proc];
      dolfin_assert(coord.size() == 18*rfacet.size());

      for (unsigned int j = 0; j < rfacet.size(); ++j)
      {
        // FIXME: inefficient? but difficult to use BBT with primitives
        // so create a small Mesh for each received prism
        Mesh prism_mesh(MPI_COMM_SELF);
        MeshEditor m_ed;
        m_ed.open(prism_mesh, mesh.topology().dim(), mesh.geometry().dim());
        m_ed.init_cells(3);
        m_ed.add_cell(0, 0,1,2,3);
        m_ed.add_cell(0, 1,2,3,4);
        m_ed.add_cell(0, 2,3,4,5);
        m_ed.init_vertices(6);
        for (unsigned int v = 0; v < 6; ++v)
          m_ed.add_vertex(v, Point(3, coord.data() + j*18 + v*3));
        m_ed.close();

        // Check all local master facets against received slave data
        for (unsigned int i = 0; i < master_facets.size(); ++i)
        {
          if (check_tet_set_collision(master_mesh, i*3, prism_mesh, 0))
          {
            const std::size_t mf = master_facets[i];
            _master_to_slave[mf].push_back(proc);
            _master_to_slave[mf].push_back(rfacet[j]);
          }
        }
      }
    }

  }
}
