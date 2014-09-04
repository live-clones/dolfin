// Copyright (C) 2005-2008 Anders Logg
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
// Modified by Garth N. Wells, 2007.
// Modified by Nuno Lopes, 2008.
// Modified by Chris Richardson, 2014.
//
// First added:  2005-12-02

#include <dolfin/common/constants.h>
#include <dolfin/common/MPI.h>
#include <dolfin/common/Timer.h>
#include <dolfin/mesh/LocalMeshData.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/mesh/MeshEditor.h>
#include "BoxMesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
BoxMesh::BoxMesh(double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 std::size_t nx, std::size_t ny, std::size_t nz)
  : Mesh(MPI_COMM_WORLD)
{
  build_distributed(x0, y0, z0, x1, y1, z1, nx, ny, nz);
}
//-----------------------------------------------------------------------------
BoxMesh::BoxMesh(MPI_Comm comm,
                 double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 std::size_t nx, std::size_t ny, std::size_t nz)
  : Mesh(comm)
{
  build_distributed(x0, y0, z0, x1, y1, z1, nx, ny, nz);
}
//-----------------------------------------------------------------------------
void BoxMesh::build_distributed(double x0, double y0, double z0,
                    double x1, double y1, double z1,
                    std::size_t nx, std::size_t ny, std::size_t nz)
{
  Timer timer("Generate Box mesh (distributed)");

  const double a = x0;
  const double b = x1;
  const double c = y0;
  const double d = y1;
  const double e = z0;
  const double f = z1;

  if (std::abs(x0 - x1) < DOLFIN_EPS || std::abs(y0 - y1) < DOLFIN_EPS
      || std::abs(z0 - z1) < DOLFIN_EPS )
  {
    dolfin_error("BoxMesh.cpp",
                 "create box",
                 "Box seems to have zero width, height or depth. Consider checking your dimensions");
  }

  if ( nx < 1 || ny < 1 || nz < 1 )
  {
    dolfin_error("BoxMesh.cpp",
                 "create box",
                 "BoxMesh has non-positive number of vertices in some dimension: number of vertices must be at least 1 in each dimension");
  }

  LocalMeshData mesh_data(this->mpi_comm());
  mesh_data.gdim = 3;
  mesh_data.tdim = 3;

  mesh_data.num_global_vertices = (nx + 1)*(ny + 1)*(nz + 1);
  const std::pair<std::size_t, std::size_t> local_vertex_range = MPI::local_range(this->mpi_comm(), 
                                                                                  mesh_data.num_global_vertices);
  const std::size_t num_local_vertices = local_vertex_range.second - local_vertex_range.first;
  mesh_data.vertex_coordinates.resize(boost::extents[num_local_vertices][mesh_data.gdim]);

  // Create vertices
  std::size_t v = 0;
  for (std::size_t vertex = local_vertex_range.first; 
       vertex != local_vertex_range.second; ++vertex)
  {
    const std::size_t ix = vertex/((ny + 1)*(nz + 1));
    const std::size_t iy = (vertex - ix*(ny + 1)*(nz + 1))/(nz + 1);
    const std::size_t iz = vertex - (ix*(ny + 1) + iy)*(nz + 1);
    
    mesh_data.vertex_coordinates[v][2] = e + (static_cast<double>(iz))*(f-e) / static_cast<double>(nz);
    mesh_data.vertex_coordinates[v][1] = c + (static_cast<double>(iy))*(d-c) / static_cast<double>(ny);
    mesh_data.vertex_coordinates[v][0] = a + (static_cast<double>(ix))*(b-a) / static_cast<double>(nx);
    ++v;
  }

  // Create tetrahedra
  const std::size_t num_cubes = nx*ny*nz;
  const std::pair<std::size_t, std::size_t> cube_range = MPI::local_range(this->mpi_comm(), num_cubes);
  mesh_data.num_global_cells = 6*num_cubes;
  const std::size_t num_local_cells = 6*(cube_range.second - cube_range.first);
  boost::multi_array<std::size_t, 2>& cells = mesh_data.cell_vertices;
  mesh_data.num_vertices_per_cell = mesh_data.tdim + 1;
  cells.resize(boost::extents[num_local_cells][mesh_data.num_vertices_per_cell]);
  mesh_data.global_cell_indices.resize(num_local_cells);

  std::size_t ci = 0;
  for (std::size_t cube = cube_range.first; cube != cube_range.second; ++cube)
  {
    const std::size_t iz = cube/(ny*nx);
    const std::size_t iy = (cube - iz*ny*nx)/nx;
    const std::size_t ix = cube - (iz*ny + iy)*nx;

    const std::size_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
    const std::size_t v1 = v0 + 1;
    const std::size_t v2 = v0 + (nx + 1);
    const std::size_t v3 = v1 + (nx + 1);
    const std::size_t v4 = v0 + (nx + 1)*(ny + 1);
    const std::size_t v5 = v1 + (nx + 1)*(ny + 1);
    const std::size_t v6 = v2 + (nx + 1)*(ny + 1);
    const std::size_t v7 = v3 + (nx + 1)*(ny + 1);

    // Note that v0 < v1 < v2 < v3 < vmid.
    cells[ci][0] = v0; cells[ci][1] = v1; cells[ci][2] = v3; cells[ci][3] = v7;
    ++ci;
    cells[ci][0] = v0; cells[ci][1] = v1; cells[ci][2] = v7; cells[ci][3] = v5;
    ++ci;
    cells[ci][0] = v0; cells[ci][1] = v5; cells[ci][2] = v7; cells[ci][3] = v4;
    ++ci;
    cells[ci][0] = v0; cells[ci][1] = v3; cells[ci][2] = v2; cells[ci][3] = v7;
    ++ci;
    cells[ci][0] = v0; cells[ci][1] = v6; cells[ci][2] = v4; cells[ci][3] = v7;
    ++ci;
    cells[ci][0] = v0; cells[ci][1] = v2; cells[ci][2] = v6; cells[ci][3] = v7;
    ++ci;
  }

  for (unsigned int i = 0; i != num_local_cells; ++i)
    mesh_data.global_cell_indices[i] = i + cube_range.first * 6;

  dolfin_assert(ci == num_local_cells);

  MeshPartitioning::build_distributed_mesh(*this, mesh_data);

  rename("mesh", "Mesh of the cuboid (a,b) x (c,d) x (e,f)");

}
//-----------------------------------------------------------------------------
void BoxMesh::build(double x0, double y0, double z0,
                    double x1, double y1, double z1,
                    std::size_t nx, std::size_t ny, std::size_t nz)
{
  Timer timer("Generate Box mesh");

  // Receive mesh according to parallel policy
  if (MPI::is_receiver(this->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*this);
    return;
  }

  const double a = x0;
  const double b = x1;
  const double c = y0;
  const double d = y1;
  const double e = z0;
  const double f = z1;

  if (std::abs(x0 - x1) < DOLFIN_EPS || std::abs(y0 - y1) < DOLFIN_EPS
      || std::abs(z0 - z1) < DOLFIN_EPS )
  {
    dolfin_error("BoxMesh.cpp",
                 "create box",
                 "Box seems to have zero width, height or depth. Consider checking your dimensions");
  }

  if ( nx < 1 || ny < 1 || nz < 1 )
  {
    dolfin_error("BoxMesh.cpp",
                 "create box",
                 "BoxMesh has non-positive number of vertices in some dimension: number of vertices must be at least 1 in each dimension");
  }

  rename("mesh", "Mesh of the cuboid (a,b) x (c,d) x (e,f)");

  // Open mesh for editing
  MeshEditor editor;
  editor.open(*this, CellType::tetrahedron, 3, 3);

  // Storage for vertex coordinates
  std::vector<double> x(3);

  // Create vertices
  editor.init_vertices_global((nx + 1)*(ny + 1)*(nz + 1), (nx + 1)*(ny + 1)*(nz + 1));
  std::size_t vertex = 0;
  for (std::size_t iz = 0; iz <= nz; iz++)
  {
    x[2] = e + (static_cast<double>(iz))*(f-e) / static_cast<double>(nz);
    for (std::size_t iy = 0; iy <= ny; iy++)
    {
      x[1] = c + (static_cast<double>(iy))*(d-c) / static_cast<double>(ny);
      for (std::size_t ix = 0; ix <= nx; ix++)
      {
        x[0] = a + (static_cast<double>(ix))*(b-a) / static_cast<double>(nx);
        editor.add_vertex(vertex, x);
        vertex++;
      }
    }
  }

  // Create tetrahedra
  editor.init_cells_global(6*nx*ny*nz, 6*nx*ny*nz);
  std::size_t cell = 0;
  std::vector<std::vector<std::size_t> > cells(6, std::vector<std::size_t>(4));
  for (std::size_t iz = 0; iz < nz; iz++)
  {
    for (std::size_t iy = 0; iy < ny; iy++)
    {
      for (std::size_t ix = 0; ix < nx; ix++)
      {
        const std::size_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
        const std::size_t v1 = v0 + 1;
        const std::size_t v2 = v0 + (nx + 1);
        const std::size_t v3 = v1 + (nx + 1);
        const std::size_t v4 = v0 + (nx + 1)*(ny + 1);
        const std::size_t v5 = v1 + (nx + 1)*(ny + 1);
        const std::size_t v6 = v2 + (nx + 1)*(ny + 1);
        const std::size_t v7 = v3 + (nx + 1)*(ny + 1);

        // Note that v0 < v1 < v2 < v3 < vmid.
        cells[0][0] = v0; cells[0][1] = v1; cells[0][2] = v3; cells[0][3] = v7;
        cells[1][0] = v0; cells[1][1] = v1; cells[1][2] = v7; cells[1][3] = v5;
        cells[2][0] = v0; cells[2][1] = v5; cells[2][2] = v7; cells[2][3] = v4;
        cells[3][0] = v0; cells[3][1] = v3; cells[3][2] = v2; cells[3][3] = v7;
        cells[4][0] = v0; cells[4][1] = v6; cells[4][2] = v4; cells[4][3] = v7;
        cells[5][0] = v0; cells[5][1] = v2; cells[5][2] = v6; cells[5][3] = v7;

        // Add cells
        std::vector<std::vector<std::size_t> >::const_iterator _cell;
        for (_cell = cells.begin(); _cell != cells.end(); ++_cell)
          editor.add_cell(cell++, *_cell);
      }
    }
  }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(this->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*this);
    return;
  }
}
//-----------------------------------------------------------------------------
