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
#include "BoxMeshDistributed.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
BoxMeshDistributed::BoxMeshDistributed(double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 std::size_t nx, std::size_t ny, std::size_t nz)
  : Mesh(MPI_COMM_WORLD)
{
  build_distributed(x0, y0, z0, x1, y1, z1, nx, ny, nz);
}
//-----------------------------------------------------------------------------
BoxMeshDistributed::BoxMeshDistributed(MPI_Comm comm,
                 double x0, double y0, double z0,
                 double x1, double y1, double z1,
                 std::size_t nx, std::size_t ny, std::size_t nz)
  : Mesh(comm)
{
  build_distributed(x0, y0, z0, x1, y1, z1, nx, ny, nz);
}
//-----------------------------------------------------------------------------
void BoxMeshDistributed::build_distributed(double x0, double y0, double z0,
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
    dolfin_error("BoxMeshDistributed.cpp",
                 "create box",
                 "Box seems to have zero width, height or depth. Consider checking your dimensions");
  }

  if ( nx < 1 || ny < 1 || nz < 1 )
  {
    dolfin_error("BoxMeshDistributed.cpp",
                 "create box",
                 "BoxMeshDistributed has non-positive number of vertices in some dimension: number of vertices must be at least 1 in each dimension");
  }

  LocalMeshData mesh_data(this->mpi_comm());
  mesh_data.gdim = 3;
  mesh_data.tdim = 3;

  mesh_data.num_global_vertices = (nx + 1)*(ny + 1)*(nz + 1);
  const std::pair<std::size_t, std::size_t> local_vertex_range =
            MPI::local_range(this->mpi_comm(),
                             mesh_data.num_global_vertices);
  const std::size_t num_local_vertices = local_vertex_range.second
                                       - local_vertex_range.first;
  mesh_data.vertex_coordinates.resize(boost::extents[num_local_vertices]
                                      [mesh_data.gdim]);

  // Create vertices
  std::size_t v = 0;
  for (std::size_t vertex = local_vertex_range.first;
       vertex != local_vertex_range.second; ++vertex)
  {
    const std::size_t nx1ny1 = (nx + 1)*(ny + 1);
    const std::size_t iz = vertex / nx1ny1;
    const std::size_t iy = (vertex % nx1ny1)/(nx + 1);
    const std::size_t ix = (vertex % nx1ny1)%(nx + 1);
    const std::size_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
    dolfin_assert(v0 == vertex);

    mesh_data.vertex_coordinates[v][2]
      = e + (static_cast<double>(iz))*(f-e) / static_cast<double>(nz);
    mesh_data.vertex_coordinates[v][1]
      = c + (static_cast<double>(iy))*(d-c) / static_cast<double>(ny);
    mesh_data.vertex_coordinates[v][0]
      = a + (static_cast<double>(ix))*(b-a) / static_cast<double>(nx);
    ++v;
  }

  // Create tetrahedra
  const std::size_t num_cubes = nx*ny*nz;

  // Find a suitable x/y/z subdivision of processes for cube
  // May crash if one dimension is small
  const std::size_t mpi_size = MPI::size(this->mpi_comm());
  const std::size_t mpi_rank = MPI::rank(this->mpi_comm());
  std::size_t npx = (int)(pow(mpi_size, 1.0/3.0) + 0.5);
  while (mpi_size%npx != 0 or npx > nx)
    npx--;
  std::size_t npy = (int)(pow(mpi_size/npx, 0.5) + 0.5);
  while ((mpi_size/npx)%npy != 0 or npy > ny)
    npy--;
  const std::size_t npz = mpi_size/npx/npy;

  std::cout << "np = " << npx << " " << npy << " " << npz <<" \n";

  const std::size_t xpos = mpi_rank/(npy*npz);
  const std::size_t ypos = (mpi_rank%(npy*npz)/npz);
  const std::size_t zpos = mpi_rank - xpos*npy*npz - ypos*npz;
  std::pair<std::size_t, std::size_t> rx = local_range(nx, xpos, npx);
  std::pair<std::size_t, std::size_t> ry = local_range(ny, ypos, npy);
  std::pair<std::size_t, std::size_t> rz = local_range(nz, zpos, npz);

  mesh_data.num_global_cells = 6*num_cubes;
  const std::size_t num_local_cells = 6*(rx.second - rx.first)
    *(ry.second - ry.first)*(rz.second - rz.first);

  boost::multi_array<std::size_t, 2>& cells = mesh_data.cell_vertices;
  mesh_data.num_vertices_per_cell = mesh_data.tdim + 1;
  cells.resize(boost::extents[num_local_cells]
               [mesh_data.num_vertices_per_cell]);
  mesh_data.global_cell_indices.resize(num_local_cells);

  std::size_t ci = 0;
  for (std::size_t ix = rx.first; ix != rx.second; ++ix)
    for (std::size_t iy = ry.first; iy != ry.second; ++iy)
      for (std::size_t iz = rz.first; iz != rz.second; ++iz)
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
	  cells[ci][0] = v0; cells[ci][1] = v1;
          cells[ci][2] = v3; cells[ci][3] = v7;
	  ++ci;
	  cells[ci][0] = v0; cells[ci][1] = v1;
          cells[ci][2] = v7; cells[ci][3] = v5;
	  ++ci;
	  cells[ci][0] = v0; cells[ci][1] = v5;
          cells[ci][2] = v7; cells[ci][3] = v4;
	  ++ci;
	  cells[ci][0] = v0; cells[ci][1] = v3;
          cells[ci][2] = v2; cells[ci][3] = v7;
	  ++ci;
	  cells[ci][0] = v0; cells[ci][1] = v6;
          cells[ci][2] = v4; cells[ci][3] = v7;
	  ++ci;
	  cells[ci][0] = v0; cells[ci][1] = v2;
          cells[ci][2] = v6; cells[ci][3] = v7;
	  ++ci;
	}

  const std::size_t cell_offset
    = MPI::global_offset(this->mpi_comm(), ci, true);

  dolfin_assert(*std::max_element(cells.data(),
				  cells.data()
				  + cells.shape()[0]*cells.shape()[1])
		< mesh_data.num_global_vertices);

  for (unsigned int i = 0; i != num_local_cells; ++i)
  {
    mesh_data.global_cell_indices[i] = i + cell_offset;
    mesh_data.cell_partition.push_back(mpi_rank);
  }

  dolfin_assert(ci == num_local_cells);

  MeshPartitioning::build_distributed_mesh(*this, mesh_data);

  rename("mesh", "Mesh of the cuboid (a,b) x (c,d) x (e,f)");
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, std::size_t>
  BoxMeshDistributed::local_range(const std::size_t N,
                                  const std::size_t block,
                                  const std::size_t nblocks) const
{
  std::pair<std::size_t, std::size_t> range;

  const std::size_t n = N / nblocks;
  const std::size_t r = N % nblocks;

  if (block < r)
    range = std::make_pair<std::size_t, std::size_t>
      (block*(n + 1), (block + 1)*(n + 1));
  else
    range = std::make_pair<std::size_t, std::size_t>
      (block*n + r, (block + 1)*n + r);

    return range;
}
//-----------------------------------------------------------------------------
