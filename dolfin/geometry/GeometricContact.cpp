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
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/common/Timer.h>

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
  // FIXME: if gdim == 2 then does Point have it's z component set to 0?
  Array<double> uval(u.value_size());

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
void GeometricContact::create_displacement_volume_mesh(Mesh& displacement_mesh,
                                                       const Mesh& mesh,
                                                       const std::vector<std::size_t> contact_facets,
                                                       const Function& u)
{
  // Precalculate triangles making the surface of a prism
  const std::size_t triangles[8][3] = {{0, 1, 2},
                                       {0, 1, 3},
                                       {1, 4, 3},
                                       {1, 2, 4},
                                       {2, 5, 4},
                                       {2, 0, 5},
                                       {0, 3, 5},
                                       {3, 4, 5}};

  // Edges of surface of prism in 2D
  const std::size_t edges[4][2] = {{0, 1},
                                       {1, 2},
                                       {2, 3},
                                       {3, 0}};

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
      for (std::size_t i = 0; i < c_per_f; ++i)
        mesh_ed.add_cell(c + i, v + triangles[i][0],
                         v + triangles[i][1],
                         v + triangles[i][2]);
    }
    else
    {
      // Add four edges
      for (std::size_t i = 0; i < c_per_f; ++i)
        mesh_ed.add_cell(c + i, v + edges[i][0],
                         v + edges[i][1]);
    }
    c += c_per_f;
    for (std::size_t i = 0; i < v_per_f; ++i)
      mesh_ed.add_vertex(v + i, master_point_set[i]);
    v += v_per_f;
  }
  mesh_ed.close();
}
//-----------------------------------------------------------------------------
void GeometricContact::create_communicated_prism_mesh(Mesh& prism_mesh,
                                                      const Mesh& mesh,
                                                      const std::vector<std::size_t>& recv_facets,
                                                      const std::vector<double>& coord)
{
  // Precalculate triangles making the surface of a prism
  const std::size_t triangles[8][3] = {{0, 1, 2},
                                       {0, 1, 3},
                                       {1, 4, 3},
                                       {1, 2, 4},
                                       {2, 5, 4},
                                       {2, 0, 5},
                                       {0, 3, 5},
                                       {3, 4, 5}};

  // Edges of surface of prism in 2D
  const std::size_t edges[4][2] = {{0, 1},
                                   {1, 2},
                                   {2, 3},
                                   {3, 0}};

  const std::size_t tdim = mesh.topology().dim();
  const std::size_t gdim = mesh.geometry().dim();

  // Find number of cells/vertices in projected prism in 2D or 3D
  const std::size_t c_per_f = GeometricContact::cells_per_facet(tdim);
  const std::size_t v_per_f = GeometricContact::vertices_per_facet(tdim);

  MeshEditor m_ed;
  m_ed.open(prism_mesh, tdim - 1, gdim);
  m_ed.init_cells(c_per_f*recv_facets.size());

  if (tdim == 3)
  {
    for (std::size_t j = 0; j < recv_facets.size(); ++j)
      for (std::size_t i = 0; i < c_per_f; ++i)
        m_ed.add_cell(j * c_per_f + i,
                      j * v_per_f + triangles[i][0],
                      j * v_per_f + triangles[i][1],
                      j * v_per_f + triangles[i][2]);
  }
  else
  {
    for (std::size_t j = 0; j < recv_facets.size(); ++j)
      for (std::size_t i = 0; i < c_per_f; ++i)
        m_ed.add_cell(j * c_per_f + i,
                      j * v_per_f + edges[i][0],
                      j * v_per_f + edges[i][1]);
  }

  m_ed.init_vertices(coord.size()/gdim);
  for (std::size_t vert = 0; vert < coord.size()/gdim; ++vert)
    m_ed.add_vertex(vert, Point(gdim, coord.data() + vert*gdim));
  m_ed.close();
}
//-----------------------------------------------------------------------------
void GeometricContact::tabulate_off_process_displacement_volume_mesh_pairs(
    const Mesh& mesh, const Mesh& slave_mesh, const Mesh& master_mesh,
    const std::vector<std::size_t>& slave_facets,
    const std::vector<std::size_t>& master_facets,
    std::map<std::size_t, std::vector<std::size_t>>& contact_facet_map)
{
  Timer ta("YYYY GC:: total inner func");
  std::vector<std::vector<std::size_t>> send_facets;
  std::vector<std::vector<double>> send_coordinates;

  const std::size_t mpi_size = MPI::size(mesh.mpi_comm());
  const std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());

  const std::size_t tdim = mesh.topology().dim();
  const std::size_t gdim = mesh.geometry().dim();

  const std::size_t c_per_f = GeometricContact::cells_per_facet(tdim);
  const std::size_t v_per_f = GeometricContact::vertices_per_facet(tdim);

  auto slave_bb = slave_mesh.bounding_box_tree();
  auto master_bb = master_mesh.bounding_box_tree();

  // Find which master processes collide with which local slave cells
  Timer t0("YYYY GC:: bbox proc colls");
  auto overlap = master_bb->compute_process_entity_collisions(*slave_bb);
  t0.stop();
  auto master_procs = overlap.first;
  auto slave_cells = overlap.second;

  // Get slave facet indices to send
  // The slave mesh consists of repeated units of eight triangles
  // making the prisms which are projected forward. The "prism index"
  // can be obtained by integer division of the cell index.
  // Each prism corresponds to a slave facet of the original mesh
  Timer t1("YYYY GC:: tabulation");
  send_facets.resize(mpi_size);
  for (std::size_t i = 0; i < master_procs.size(); ++i)
  {
    const std::size_t master_rank = master_procs[i];
    // Ignore local (already done)
    if (master_rank != mpi_rank)
    {
      // Get facet from cell index (2D: 4 edges per prism, 3D: eight triangles per prism)
      std::size_t facet = slave_cells[i]/c_per_f;
      send_facets[master_rank].push_back(facet);
    }
  }

  // Get unique set of facets to send to each process
  for (std::size_t p = 0; p != mpi_size; ++p)
  {
    std::vector<std::size_t>& v = send_facets[p];
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
  }

  // Get coordinates of prism (18 doubles in 3D, 8 in 2D)
  // and convert index of slave mesh back to index of main mesh
  send_coordinates.resize(mpi_size);
  for (std::size_t p = 0; p != mpi_size; ++p)
  {
    std::vector<double>& coords_p = send_coordinates[p];
    for (auto& q : send_facets[p])
    {
      const double* coord_vals = slave_mesh.geometry().x(q*v_per_f);
      coords_p.insert(coords_p.end(), coord_vals, coord_vals + gdim*v_per_f);
      // Convert from local prism index to main mesh local facet indexing
      q = slave_facets[q];
    }
  }
  t1.stop();

  Timer t5("YYYY GC:: all_to_all all the things");
  std::vector<std::vector<std::size_t>> recv_facets(mpi_size);
  std::vector<std::vector<double>> recv_coordinates(mpi_size);

  MPI::all_to_all(mesh.mpi_comm(), send_facets, recv_facets);
  MPI::all_to_all(mesh.mpi_comm(), send_coordinates, recv_coordinates);
  t5.stop();

  if (master_mesh.num_vertices() == 0)
    return;

  Mesh master_mesh_on_proc(MPI_COMM_SELF);
  GeometricContact::create_on_process_sub_mesh(master_mesh_on_proc, master_mesh);

  Timer t6("YYYY GC:: Unpack everything");
  for (std::size_t proc = 0; proc != mpi_size; ++proc)
  {
    const auto& rfacet = recv_facets[proc];
    const auto& coord = recv_coordinates[proc];
    dolfin_assert(coord.size() == gdim*v_per_f*rfacet.size());

    // Received facets is empty, so there can be no collisions on this process
    if (rfacet.empty())
      continue;

    Mesh prism_mesh(MPI_COMM_SELF);
    GeometricContact::create_communicated_prism_mesh(prism_mesh, mesh, rfacet, coord);
    GeometricContact::tabulate_on_process_bbox_collisions(proc, master_mesh_on_proc, master_facets,
                                                          prism_mesh, rfacet, contact_facet_map);
  }
  t6.stop();
  ta.stop();
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

  Timer t("XXXX GC:: create_displacement_volume_mesh");
  // Make the displacement volume mesh of the master facets
  Mesh master_mesh(mesh.mpi_comm());
  GeometricContact::create_displacement_volume_mesh(master_mesh, mesh, master_facets, u);

  // Make the displacement volume mesh of the slave facets
  Mesh slave_mesh(mesh.mpi_comm());
  GeometricContact::create_displacement_volume_mesh(slave_mesh, mesh, slave_facets, u);
  t.stop();

  _master_to_slave.clear();
  _slave_to_master.clear();

  std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());
  std::size_t mpi_size = MPI::size(mesh.mpi_comm());

  // Find number of cells/vertices in projected prism in 2D or 3D
  const std::size_t c_per_f = GeometricContact::cells_per_facet(tdim);

  // Check each master 'prism' against each slave 'prism'
  // Map is stored as local_master_facet -> [mpi_rank, local_index, mpi_rank, local_index ...]
  // First check locally
  Timer t1("XXXX GC:: collision detec");
  if (master_mesh.num_vertices() != 0 and slave_mesh.num_vertices() != 0)
  {
    Mesh master_mesh_on_proc(MPI_COMM_SELF);
    GeometricContact::create_on_process_sub_mesh(master_mesh_on_proc, master_mesh);
    Mesh slave_mesh_on_proc(MPI_COMM_SELF);
    GeometricContact::create_on_process_sub_mesh(slave_mesh_on_proc, slave_mesh);
    GeometricContact::tabulate_on_process_bbox_collisions(mpi_rank,
                                                          master_mesh_on_proc, master_facets,
                                                          slave_mesh_on_proc, slave_facets,
                                                          _master_to_slave);
    GeometricContact::tabulate_on_process_bbox_collisions(mpi_rank,
                                                          slave_mesh_on_proc, slave_facets,
                                                          master_mesh_on_proc, master_facets,
                                                          _slave_to_master);
  }
  t1.stop();

  // Find which [master global/slave entity] BBs overlap in parallel
  if (mpi_size > 1)
  {
    Timer t("XXXX GC:: tabulate_off_process_displacement_volume_mesh_pairs");
    GeometricContact::tabulate_off_process_displacement_volume_mesh_pairs(mesh, slave_mesh, master_mesh, slave_facets,
                                                                          master_facets, _master_to_slave);

    GeometricContact::tabulate_off_process_displacement_volume_mesh_pairs(mesh, master_mesh, slave_mesh, master_facets,
                                                                          slave_facets, _slave_to_master);
    t.stop();
  }
}
//-----------------------------------------------------------------------------
void GeometricContact::tabulate_collided_cell_dofs(const Mesh& mesh, const GenericDofMap& dofmap,
                                                   const std::map<std::size_t, std::vector<std::size_t>>& master_to_slave,
                                                   std::map<std::size_t, std::vector<std::size_t>>& facet_to_contacted_dofs,
                                                   std::map<std::size_t, std::vector<std::size_t>>& facet_to_off_proc_contacted_dofs)
{
  const std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());
  const std::size_t mpi_size = MPI::size(mesh.mpi_comm());

  // tabulate global DoFs required for off-process insertion
  std::vector<std::size_t> local_to_global_dofs;
  dofmap.tabulate_local_to_global_dofs(local_to_global_dofs);

  const std::size_t tdim = mesh.topology().dim();

  // Send the master cell's dofs to the slave.
  // [proc: [local_slave, contact master dofs, local slave, contact master dofs, ...]]
  std::vector<std::vector<std::size_t>> send_master_dofs(mpi_size);

  Timer t0("VVVVV GC: tabulate dofs to send");
  for (const auto& m2s : master_to_slave)
  {
    const std::size_t mi = m2s.first;  // master facet indices
    const std::vector<std::size_t> sis = m2s.second;  // [proc, slave facet idx, proc, slave facet idx, ...]

    // Cell to which the master facet belongs and its DoFs
    const Cell m_cell(mesh, Facet(mesh, mi).entities(tdim)[0]);
    const auto& m_cell_dofs_eigen_map
        = dofmap.cell_dofs(m_cell.index());
    const auto m_cell_dofs = ArrayView<const dolfin::la_index>(
        m_cell_dofs_eigen_map.size(), m_cell_dofs_eigen_map.data());

    for (std::size_t j=0; j<sis.size(); j+=2)
    {
      const std::size_t& slave_proc = sis[j];
      const std::size_t& si = sis[j+1];

      // If the slave is on the current process, then append the dofs.
      if (mpi_rank == slave_proc)
      {
        // Cell to which the slave facet belongs and its DoFs
        const Cell s_cell(mesh, Facet(mesh, si).entities(tdim)[0]);
        const auto& s_cell_dofs_eigen_map
            = dofmap.cell_dofs(s_cell.index());
        const auto s_cell_dofs = ArrayView<const dolfin::la_index>(
            s_cell_dofs_eigen_map.size(), s_cell_dofs_eigen_map.data());

        // Insert slave dofs into the map
//        facet_to_contacted_dofs[mi].insert(std::end(facet_to_contacted_dofs[mi]),
//                                           std::begin(s_cell_dofs),
//                                           std::end(s_cell_dofs));

        facet_to_contacted_dofs[si].insert(std::end(facet_to_contacted_dofs[si]),
                                           std::begin(m_cell_dofs),
                                           std::end(m_cell_dofs));
      }
      else // schedule master dofs to be dispatched to the slaves
      {
        // Convert to global dofs for column dof entry
        std::vector<std::size_t> global_master_dofs(std::begin(m_cell_dofs),
                                                    std::end(m_cell_dofs));
        for (auto& dof : global_master_dofs)
          dof = local_to_global_dofs[dof];

        send_master_dofs[slave_proc].push_back(si);
        send_master_dofs[slave_proc].insert(std::end(send_master_dofs[slave_proc]),
                                            std::begin(global_master_dofs),
                                            std::end(global_master_dofs));
      }
    }

  }
  t0.stop();

  // Running in serial. So we're done.
  if (mpi_size == 1)
    return;

  Timer t1("VVVVV GC: tabulate dofs all to all");
  std::vector<std::size_t> recv_master_dofs;
  MPI::all_to_all(mesh.mpi_comm(), send_master_dofs, recv_master_dofs);
  t1.stop();

  const std::size_t num_dofs_per_cell = dofmap.max_element_dofs();

  // Tabulate the communicated dofs belonging to the master cells. This requires
  // saving the global DoFs of the master cells to _local_facet_to_off_proc_contact_dofs
  // at the index of the local on process slave cell index.
  Timer t3("VVVVV GC: get dofs out");
  for (std::size_t j = 0; j < recv_master_dofs.size(); j += num_dofs_per_cell + 1)
  {
    // global_master_dofs[j] is the slave facet local index. global_master_dofs[j+1:j+num_dofs_per_cell]
    // are the communicated master cell dofs.
    std::vector<std::size_t>& dof_list = facet_to_off_proc_contacted_dofs[recv_master_dofs[j]];
    for (std::size_t i = 0; i < num_dofs_per_cell; ++i)
      dof_list.push_back(recv_master_dofs[j + i + 1]);
  }
  t3.stop();
}
//-----------------------------------------------------------------------------
void GeometricContact::tabulate_contact_cell_to_shared_dofs(const Mesh& mesh, const Function& u,
                                                            const std::vector<std::size_t>& master_facets,
                                                            const std::vector<std::size_t>& slave_facets)
{
  const auto V = u.function_space();
  const auto dofmap = V->dofmap();

  // Start from fresh
  _local_facet_to_contact_dofs.clear();
  _local_facet_to_off_proc_contact_dofs.clear();

  GeometricContact::tabulate_collided_cell_dofs(mesh, *dofmap, _master_to_slave,
                                                _local_facet_to_contact_dofs,
                                                _local_facet_to_off_proc_contact_dofs);
  GeometricContact::tabulate_collided_cell_dofs(mesh, *dofmap, _slave_to_master,
                                                _local_facet_to_contact_dofs,
                                                _local_facet_to_off_proc_contact_dofs);
}
//-----------------------------------------------------------------------------
void
GeometricContact::tabulate_contact_shared_cells(Mesh& mesh, Function& u,
                                                const std::vector<std::size_t>& master_facets,
                                                const std::vector<std::size_t>& slave_facets)
{
  const std::size_t mpi_size = MPI::size(mesh.mpi_comm());

  const auto& s2m = slave_to_master();

  const auto V = u.function_space();
  const auto dofmap = V->dofmap();

  // tabulate global DoFs required for off-process insertion
  std::vector<std::size_t> local_to_global_dofs;
  dofmap->tabulate_local_to_global_dofs(local_to_global_dofs);

  const std::size_t tdim = mesh.topology().dim();

  // Communicate the slave cells' meta data to the master.
  // First we tabulate the slave cells' information for dispatch
  std::vector<std::vector<std::size_t>> slave_facet_infos_send(mpi_size);
  std::vector<std::vector<std::size_t>> slave_cell_global_dofs_send(mpi_size);
  std::vector<std::vector<double>> slave_cell_dof_coords_send(mpi_size);
  std::vector<std::vector<double>> slave_cell_dof_coeffs_send(mpi_size);

  std::vector<double> dof_coords;

  const FiniteElement element = *V->element();
  std::vector<double> dof_coeffs(element.space_dimension());
  ufc::cell ufc_cell;

  for (const auto& slave : s2m)
  {
    const std::size_t slave_idx = slave.first; // slave facet indices
    const std::vector<std::size_t> master_procs_idxs = slave.second;  // [proc, master facet idx, proc, master facet idx, ...

    const Facet slave_facet(mesh, slave_idx);
    const Cell slave_cell = Cell(mesh, slave_facet.entities(tdim)[0]);
    const std::size_t slave_facet_local_idx = slave_cell.index(slave_facet);

    // Tabulate dof coords
    slave_cell.get_coordinate_dofs(dof_coords);

    // Tabulate cell dofs
    const auto& s_cell_dofs_eigen_map
        = dofmap->cell_dofs(slave_cell.index());
    const auto s_cell_dofs = ArrayView<const dolfin::la_index>(
        s_cell_dofs_eigen_map.size(), s_cell_dofs_eigen_map.data());
    // Convert to global dofs
    std::vector<std::size_t> global_s_cell_dofs(std::begin(s_cell_dofs),
                                                std::end(s_cell_dofs));
    for (auto& dof : global_s_cell_dofs)
      dof = local_to_global_dofs[dof];

    // Tabulate cell dof coefficients (local finite element vector on the slave element)

    // Update to current cell
    slave_cell.get_cell_data(ufc_cell);
    u.restrict(dof_coeffs.data(), element, slave_cell, dof_coords.data(),
             ufc_cell);

    for (std::size_t j=0; j<master_procs_idxs.size(); j+=2)
    {
      const std::size_t master_proc = master_procs_idxs[j];
      const std::size_t master_facet_idx = master_procs_idxs[j+1];

      // [slave_facet_idx (local index to mesh), slave_facet_local_idx (local to the cell), master_facet_idx]
      slave_facet_infos_send[master_proc].push_back(slave_facet.index());
      slave_facet_infos_send[master_proc].push_back(slave_facet_local_idx);
      slave_facet_infos_send[master_proc].push_back(master_facet_idx);

      slave_cell_global_dofs_send[master_proc].insert(std::end(slave_cell_global_dofs_send[master_proc]),
                                                      std::begin(global_s_cell_dofs), std::end(global_s_cell_dofs));

      slave_cell_dof_coords_send[master_proc].insert(std::end(slave_cell_dof_coords_send[master_proc]),
                                                     std::begin(dof_coords), std::end(dof_coords));

      slave_cell_dof_coeffs_send[master_proc].insert(std::end(slave_cell_dof_coeffs_send[master_proc]),
                                                      std::begin(dof_coeffs), std::end(dof_coeffs));
    }
  }

  // Receive the cell metadata
  std::vector<std::vector<std::size_t>> slave_facet_infos_recv;
  std::vector<std::vector<std::size_t>> slave_cell_global_dofs_recv;
  std::vector<std::vector<double>> slave_cell_dof_coords_recv;
  std::vector<std::vector<double>> slave_cell_dof_coeffs_recv;

  // FIXME: This is too much MP....?
  MPI::all_to_all(mesh.mpi_comm(), slave_facet_infos_send, slave_facet_infos_recv);
  MPI::all_to_all(mesh.mpi_comm(), slave_cell_global_dofs_send, slave_cell_global_dofs_recv);
  MPI::all_to_all(mesh.mpi_comm(), slave_cell_dof_coords_send, slave_cell_dof_coords_recv);
  MPI::all_to_all(mesh.mpi_comm(), slave_cell_dof_coeffs_send, slave_cell_dof_coeffs_recv);

  const std::size_t num_coords_per_cell = Cell(mesh, 0).num_vertices()*mesh.geometry().dim();
  const std::size_t num_dofs_per_cell = dofmap->max_element_dofs();

  for (std::size_t proc_source=0; proc_source<mpi_size; ++proc_source)
  {
    if (slave_facet_infos_recv[proc_source].empty())
      continue;

    for (std::size_t j=0; j<slave_facet_infos_recv[proc_source].size()/3; ++j)
    {
      // [slave_facet_idx, slave_facet_local_idx, master_facet_idx]
      const std::size_t slave_facet_idx = slave_facet_infos_recv[proc_source][3*j];
      const std::size_t slave_facet_local_idx = slave_facet_infos_recv[proc_source][3*j + 1];
      const std::size_t master_idx = slave_facet_infos_recv[proc_source][3*j + 2];

      const std::vector<std::size_t> cell_dofs(slave_cell_global_dofs_recv[proc_source].begin() + j*num_dofs_per_cell,
                                               slave_cell_global_dofs_recv[proc_source].begin() + (j + 1)*num_dofs_per_cell);
      const std::vector<double> cell_dof_coords(slave_cell_dof_coords_recv[proc_source].begin() + j*num_coords_per_cell,
                                                slave_cell_dof_coords_recv[proc_source].begin() + (j + 1)*num_coords_per_cell);

      // FIXME: Contact on manifolds will not have the same number of coeffs as coords
      const std::vector<double> cell_dof_coeffs(slave_cell_dof_coeffs_recv[proc_source].begin() + j*num_coords_per_cell,
                                                slave_cell_dof_coeffs_recv[proc_source].begin() + (j + 1)*num_coords_per_cell);

      auto cell_md = std::make_shared<CellMetaData>(slave_facet_idx, slave_facet_local_idx,
                           cell_dof_coords, cell_dofs, cell_dof_coeffs);
      _master_facet_to_contacted_cells[master_idx].push_back(cell_md);
    }
  }

}
//-----------------------------------------------------------------------------
void GeometricContact::tabulate_on_process_bbox_collisions(const std::size_t mpi_rank,
                                                           const Mesh& master_mesh,
                                                           const std::vector<std::size_t>& master_facets,
                                                           const Mesh& slave_mesh,
                                                           const std::vector<std::size_t>& slave_facets,
                                                           std::map<std::size_t, std::vector<std::size_t>>& master_to_slave)
{
  // Master and slave meshes are tdim = D-1 manifolds in gdim = D.
  const std::size_t c_per_f = GeometricContact::cells_per_facet(master_mesh.topology().dim() + 1);

  BoundingBoxTree bbox_master;
  bbox_master.build(master_mesh);

  BoundingBoxTree bbox_slave;
  bbox_slave.build(slave_mesh);

  const auto master_slave_bbox_coll =
      bbox_master.compute_entity_collisions(bbox_slave);

  const auto& master_idxs = master_slave_bbox_coll.first;
  const auto& slave_idxs = master_slave_bbox_coll.second;

  std::map<std::size_t, std::set<std::size_t>> m2s_set;
  for (std::size_t j=0; j<master_idxs.size(); ++j)
  {
    const std::size_t mf = master_facets[master_idxs[j]/c_per_f];
    const std::size_t sf = slave_facets[slave_idxs[j]/c_per_f];
    m2s_set[mf].insert(sf);
  }

  for(const auto& m2s_pair : m2s_set)
  {
    const std::size_t mf = m2s_pair.first;
    const std::set<std::size_t>& sfs = m2s_pair.second;
    for (const auto& sf : sfs)
    {
      master_to_slave[mf].push_back(mpi_rank);
      master_to_slave[mf].push_back(sf);
    }
  }
}
//-----------------------------------------------------------------------------
void GeometricContact::create_on_process_sub_mesh(Mesh& sub_mesh, const Mesh& mesh)
{
  MeshEditor me;
  me.open(sub_mesh, mesh.topology().dim(), mesh.geometry().dim());

  me.init_cells(mesh.num_cells());
  me.init_vertices(mesh.num_vertices());

  const std::size_t v_per_c = mesh.topology().dim() + 1;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    const auto v_idxs_ptr = c->entities(0);
    // FIXME: The interface to me.add_cell is frustrating for variable length arguments
    // Need to fix this so we don't instantiate a std::vector each iteration.
    std::vector<std::size_t> v_idxs(v_idxs_ptr, v_idxs_ptr + v_per_c);
    me.add_cell(c->index(), v_idxs);
  }
  for (VertexIterator v(mesh); !v.end(); ++v)
    me.add_vertex(v->index(), v->point());
  me.close();
}
//-----------------------------------------------------------------------------
