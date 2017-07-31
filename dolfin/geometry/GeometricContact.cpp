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
std::vector<Point> GeometricContact::create_deformed_segment_volume(const Mesh& mesh,
                                                                    std::size_t facet_index,
                                                                    const Function& u,
                                                                    std::size_t gdim)
{
  if (gdim != 3 and gdim != 2)
  {
    dolfin_error("GeometricContact.cpp",
                 "calculate deformation volume",
                 "Must be 2D or 3D");
  }

  Facet facet(mesh, facet_index);

  // Get id of attached cell
  std::size_t id = facet.entities(mesh.topology().dim())[0];

  // Vector value of Function
  Array<double> uval(3);
  // make sure z-component is zero, for 2D
  uval[2] = 0.0;

  const Cell cell(mesh, id);
  ufc::cell ufc_cell;
  cell.get_cell_data(ufc_cell);

  Point X1 = Vertex(mesh, facet.entities(0)[0]).point();
  const Array<double> _x1(3, X1.coordinates());
  u.eval(uval, _x1, cell, ufc_cell);
  Point x1 = X1 + Point(uval);

  Point X2 = Vertex(mesh, facet.entities(0)[1]).point();
  const Array<double> _x2(3, X2.coordinates());
  u.eval(uval, _x2, cell, ufc_cell);
  Point x2 = X2 + Point(uval);

  if (gdim == 2)
    return {X1, X2, x1, x2};

  Point X3 = Vertex(mesh, facet.entities(0)[2]).point();
  const Array<double> _x3(3, X3.coordinates());
  u.eval(uval, _x3, cell, ufc_cell);
  Point x3 = X3 + Point(uval);

  return {X1, X2, X3, x1, x2, x3};
}
//-----------------------------------------------------------------------------
bool GeometricContact::check_tri_set_collision(const Mesh& master_mesh, std::size_t mindex,
                                               const Mesh& slave_mesh, std::size_t sindex)
{

  for (unsigned int i = mindex; i < mindex + 8; ++i)
    for (unsigned int j = sindex; j < sindex + 8; ++j)
    {
      if (CollisionDetection::collides_triangle_triangle(Cell(master_mesh, i),
                                                         Cell(slave_mesh, j)))
        return true;
    }

  return false;
}
//-----------------------------------------------------------------------------
bool GeometricContact::check_edge_set_collision(const Mesh& master_mesh, std::size_t mindex,
                                               const Mesh& slave_mesh, std::size_t sindex)
{
  const auto mconn = master_mesh.topology()(1, 0);
  const auto sconn = slave_mesh.topology()(1, 0);

  for (unsigned int i = mindex; i < mindex + 4; ++i)
    for (unsigned int j = sindex; j < sindex + 4; ++j)
    {
      const Point p0 = Vertex(master_mesh, mconn(i)[0]).point();
      const Point p1 = Vertex(master_mesh, mconn(i)[1]).point();
      const Point p2 = Vertex(slave_mesh, sconn(j)[0]).point();
      const Point p3 = Vertex(slave_mesh, sconn(j)[1]).point();

      if (CollisionDetection::collides_edge_edge(p0, p1, p2, p3))
        return true;
    }

  return false;
}
//-----------------------------------------------------------------------------
bool GeometricContact::create_displacement_volume_mesh(Mesh& displacement_mesh,
                                                       const Mesh& mesh,
                                                       const std::vector<std::size_t> contact_facets,
                                                       const Function& u)
{
  const std::size_t tdim = mesh.topology().dim();

  // Find number of cells/vertices in projected prism in 2D or 3D
  const std::size_t c_per_f = GeometricContact::cells_per_facet(tdim);
  const std::size_t v_per_f = GeometricContact::vertices_per_facet(tdim);

  // Local mesh of master 'prisms', eight triangles are created per facet in 3D
  // Four edges in 2D
  MeshEditor mesh_ed;
  mesh_ed.open(displacement_mesh, mesh.topology().dim() - 1, mesh.geometry().dim());
  std::size_t nf_local = contact_facets.size();
  std::size_t nf_global = MPI::sum(mesh.mpi_comm(), nf_local);
  mesh_ed.init_cells_global(nf_local*c_per_f, nf_global*c_per_f);
  mesh_ed.init_vertices_global(nf_local*v_per_f, nf_global*v_per_f);

  std::size_t c = 0;
  std::size_t v = 0;

  for (const auto& mf : contact_facets)
  {
    std::vector<Point> master_point_set = create_deformed_segment_volume(mesh, mf, u, tdim);

    if (tdim == 3)
    {
      // Add eight triangles
      for (unsigned int i = 0; i < c_per_f; ++i)
        mesh_ed.add_cell(c + i, v + triangles[i][0],
                         v + triangles[i][1],
                         v + triangles[i][2]);
    }
    else
    {
      // Add four edges
      for (unsigned int i = 0; i < c_per_f; ++i)
        mesh_ed.add_cell(c + i, v + edges[i][0],
                         v + edges[i][1]);
    }
    c += c_per_f;
    for (unsigned int i = 0; i < v_per_f; ++i)
      mesh_ed.add_vertex(v + i, master_point_set[i]);
    v += v_per_f;
  }
  mesh_ed.close();
}
//-----------------------------------------------------------------------------
bool GeometricContact::create_communicated_prism_mesh(Mesh& prism_mesh,
                                                      const Mesh& mesh,
                                                      const std::vector<std::size_t>& facet,
                                                      const std::vector<double>& coord,
                                                      std::size_t local_facet_idx)
{
  const std::size_t tdim = mesh.topology().dim();
  const std::size_t gdim = mesh.geometry().dim();

  // Find number of cells/vertices in projected prism in 2D or 3D
  const std::size_t c_per_f = GeometricContact::cells_per_facet(tdim);
  const std::size_t v_per_f = GeometricContact::vertices_per_facet(tdim);

  MeshEditor m_ed;
  m_ed.open(prism_mesh, tdim - 1, gdim);
  m_ed.init_cells(c_per_f);
  if (tdim == 3)
  {
    for (unsigned int i = 0; i < c_per_f; ++i)
      m_ed.add_cell(i, triangles[i][0], triangles[i][1], triangles[i][2]);
  }
  else
  {
    for (unsigned int i = 0; i < c_per_f; ++i)
      m_ed.add_cell(i, edges[i][0], edges[i][1]);
  }

  m_ed.init_vertices(v_per_f);
  for (unsigned int vert = 0; vert < v_per_f; ++vert)
    m_ed.add_vertex(vert, Point(gdim, coord.data() + (local_facet_idx*v_per_f + vert)*gdim));
  m_ed.close();
}
//-----------------------------------------------------------------------------
void GeometricContact::contact_surface_map_volume_sweep(Mesh& mesh, Function& u,
const std::vector<std::size_t>& master_facets, const std::vector<std::size_t>& slave_facets)
{
  // Construct a dictionary mapping master facets to their collided slave counterparts.

  const std::size_t tdim = mesh.topology().dim();
  if (tdim != 3 and tdim != 2)
  {
    dolfin_error("GeometricContact.cpp",
                 "find contact surface",
                 "Only implemented in 2D/3D");
  }

  const std::size_t gdim = mesh.geometry().dim();
  if (gdim != tdim)
  {
    dolfin_error("GeometricContact.cpp",
                 "find contact surface",
                 "Manifold meshes not supported");
  }

  // Ensure BBT built
  auto mesh_bb = mesh.bounding_box_tree();

  // Make sure facet->cell connections are made
  mesh.init(tdim - 1, tdim);

  Mesh master_mesh(mesh.mpi_comm());
  GeometricContact::create_displacement_volume_mesh(master_mesh, mesh, master_facets, u);
  auto master_bb = master_mesh.bounding_box_tree();

  Mesh slave_mesh(mesh.mpi_comm());
  GeometricContact::create_displacement_volume_mesh(slave_mesh, mesh, slave_facets, u);
  auto slave_bb = slave_mesh.bounding_box_tree();

  _master_to_slave.clear();

  std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());
  std::size_t mpi_size = MPI::size(mesh.mpi_comm());

  // Find number of cells/vertices in projected prism in 2D or 3D
  const std::size_t cells_per_facet = (tdim - 1)*4;
  const std::size_t vertices_per_facet = tdim*2;

  // Check each master 'prism' against each slave 'prism'
  // Map is stored as local_master_facet -> [mpi_rank, local_index, mpi_rank, local_index ...]
  // First check locally
  for (unsigned int i = 0; i < master_facets.size(); ++i)
    for (unsigned int j = 0; j < slave_facets.size(); ++j)
    {
      // FIXME: for efficiency, use BBT here
      bool collision;
      if (tdim == 3)
        collision = check_tri_set_collision(master_mesh, i*cells_per_facet, slave_mesh, j*cells_per_facet);
      else
        collision = check_edge_set_collision(master_mesh, i*cells_per_facet, slave_mesh, j*cells_per_facet);

      if (collision)
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
    // The slave mesh consists of repeated units of eight triangles
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
        // Get facet from cell index (eight triangles per prism)
        unsigned int facet = slave_cells[i]/cells_per_facet;
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

    // Get coordinates of prism (18 doubles in 3D, 8 in 2D)
    // and convert index of slave mesh back to index of main mesh
    std::vector<std::vector<double>> send_coordinates(mpi_size);
    for (unsigned int p = 0; p != mpi_size; ++p)
    {
      std::vector<double>& coords_p = send_coordinates[p];
      for (auto& q : send_facets[p])
      {
        const double *coord_vals = slave_mesh.geometry().vertex_coordinates(q*vertices_per_facet);
        coords_p.insert(coords_p.end(), coord_vals, coord_vals + gdim*vertices_per_facet);
        // Convert back to main mesh indexing
        q = slave_facets[q];
      }
    }

    std::vector<std::vector<std::size_t>> recv_facets(mpi_size);
    MPI::all_to_all(mesh.mpi_comm(), send_facets, recv_facets);

    std::vector<std::vector<double>> recv_coordinates(mpi_size);
    MPI::all_to_all(mesh.mpi_comm(), send_coordinates, recv_coordinates);

    for (std::size_t proc = 0; proc != mpi_size; ++proc)
    {
      auto& rfacet = recv_facets[proc];
      auto& coord = recv_coordinates[proc];
      dolfin_assert(coord.size() == gdim*vertices_per_facet*rfacet.size());

      for (unsigned int j = 0; j < rfacet.size(); ++j)
      {
        // FIXME: inefficient? but difficult to use BBT with primitives
        // so create a small Mesh for each received prism
        Mesh prism_mesh(MPI_COMM_SELF);
        GeometricContact::create_communicated_prism_mesh(prism_mesh, mesh, rfacet, coord, j);

        // Check all local master facets against received slave data
        for (unsigned int i = 0; i < master_facets.size(); ++i)
        {
          bool collision;
          if (tdim == 3)
            collision = check_tri_set_collision(master_mesh, i*cells_per_facet, prism_mesh, 0);
          else
            collision = check_edge_set_collision(master_mesh, i*cells_per_facet, prism_mesh, 0);

          if (collision)
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
////-----------------------------------------------------------------------------
//void
//tabulate_contact_cell_to_shared_dofs(Mesh& mesh, Function& u,
//                                     const std::vector<std::size_t>& master_facets,
//                                     const std::vector<std::size_t>& slave_facets)
//{
//
//}
