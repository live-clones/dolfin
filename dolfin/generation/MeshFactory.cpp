// Copyright (C) 2005-2015 Anders Logg, Chris Richardson
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
// Modified by Garth N. Wells 2007
// Modified by Nuno Lopes 2008
// Modified by Kristian B. Oelgaard 2009
// Modified by Chris Richardson 2015

#include <dolfin/common/Timer.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/mesh/Vertex.h>

#include "MeshFactory.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::UnitCubeMesh(MPI_Comm mpi_comm,
                                                std::size_t nx,
                                                std::size_t ny,
                                                std::size_t nz,
                                                std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));
  build_box_mesh(mesh, Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0),
                 nx, ny, nz, options);
  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::BoxMesh(MPI_Comm mpi_comm,
                                           const Point& p0, const Point& p1,
                                           std::size_t nx,
                                           std::size_t ny,
                                           std::size_t nz,
                                           std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));
  build_box_mesh(mesh, p0, p1, nx, ny, nz, options);
  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::UnitSquareMesh(MPI_Comm mpi_comm,
                                                  std::size_t nx,
                                                  std::size_t ny,
                                                  std::string options)
{
  return MeshFactory::RectangleMesh(mpi_comm, Point(0.0, 0.0), Point(1.0, 1.0),
                                    nx, ny, options);
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::RectangleMesh(MPI_Comm mpi_comm,
                                                 const Point& p0,
                                                 const Point& p1,
                                                 std::size_t nx,
                                                 std::size_t ny,
                                                 std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));

  if (options != "left/right" and
      options != "right/left" and
      options != "crossed" and
      options != "left" and
      options != "right")
    dolfin_error("MeshFactory.cpp",
                 "determine mesh options",
                 "Unknown mesh diagonal definition: allowed options are \"left\", \"right\", \"left/right\" and \"crossed\"");
  build_rectangle_mesh(mesh, p0, p1, nx, ny, options);
  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::IntervalMesh(MPI_Comm mpi_comm,
                                                std::size_t nx, double a,
                                                double b, std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));

  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return mesh;
  }

  if (std::abs(a - b) < DOLFIN_EPS)
  {
    dolfin_error("Interval.cpp",
                 "create interval",
                 "Length of interval is zero. Consider checking your dimensions");
  }

  if (b < a)
  {
    dolfin_error("Interval.cpp",
                 "create interval",
                 "Length of interval is negative. Consider checking the order of your arguments");
  }

  if (nx < 1)
  {
    dolfin_error("Interval.cpp",
                 "create interval",
                 "Number of points on interval is (%d), it must be at least 1", nx);
  }

  rename("mesh", "Mesh of the interval (a, b)");

  // Open mesh for editing
  MeshEditor editor;
  editor.open(*mesh, CellType::interval, 1, 1);

  // Create vertices and cells:
  editor.init_vertices_global((nx+1), (nx+1));
  editor.init_cells_global(nx, nx);

  // Create main vertices:
  for (std::size_t ix = 0; ix <= nx; ix++)
  {
    const std::vector<double>
      x(1, a + (static_cast<double>(ix)*(b - a)/static_cast<double>(nx)));
    editor.add_vertex(ix, x);
  }

  // Create intervals
  for (std::size_t ix = 0; ix < nx; ix++)
  {
    std::vector<std::size_t> cell(2);
    cell[0] = ix; cell[1] = ix + 1;
    editor.add_cell(ix, cell);
  }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh->mpi_comm()))
  {
    std::cout << "Building mesh (dist 0a)" << std::endl;
    MeshPartitioning::build_distributed_mesh(*mesh);
    std::cout << "Building mesh (dist 1a)" << std::endl;
  }

  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh>
MeshFactory::UnitQuadMesh(MPI_Comm mpi_comm, std::size_t nx, std::size_t ny,
                          std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));

  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return mesh;
  }

  MeshEditor editor;
  editor.open(*mesh, CellType::quadrilateral, 2, 2);

  // Create vertices and cells:
  editor.init_vertices_global((nx + 1)*(ny + 1), (nx + 1)*(ny + 1));
  editor.init_cells_global(nx*ny, nx*ny);

  // Storage for vertices
  std::vector<double> x(2);

  const double a = 0.0;
  const double b = 1.0;
  const double c = 0.0;
  const double d = 1.0;

  // Create main vertices:
  std::size_t vertex = 0;
  for (std::size_t iy = 0; iy <= ny; iy++)
  {
    x[1] = c + ((static_cast<double>(iy))*(d - c)/static_cast<double>(ny));
    for (std::size_t ix = 0; ix <= nx; ix++)
    {
      x[0] = a + ((static_cast<double>(ix))*(b - a)/static_cast<double>(nx));
      editor.add_vertex(vertex, x);
      vertex++;
    }
  }

  // Create rectangles
  std::size_t cell = 0;
  std::vector<std::size_t> v(4);
  for (std::size_t iy = 0; iy < ny; iy++)
    for (std::size_t ix = 0; ix < nx; ix++)
    {
      v[0] = iy*(nx + 1) + ix;
      v[1] = v[0] + 1;
      v[2] = v[0] + (nx + 1);
      v[3] = v[1] + (nx + 1);
      editor.add_cell(cell, v);
      ++cell;
    }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh->mpi_comm()))
    MeshPartitioning::build_distributed_mesh(*mesh);

  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::UnitTetrahedronMesh(MPI_Comm mpi_comm,
                                                       std::string options)
{
  if (MPI::size(mpi_comm) != 1)
    dolfin_error("MeshFactory.cpp",
                 "generate UnitTetraHedronMesh",
                 "Cannot generate distributed mesh");

  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));

  // Open mesh for editing
  MeshEditor editor;
  editor.open(*mesh, CellType::tetrahedron, 3, 3);

  // Create vertices
  editor.init_vertices_global(4, 4);
  std::vector<double> x(3);
  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
  editor.add_vertex(0, x);

  x[0] = 1.0; x[1] = 0.0; x[2] = 0.0;
  editor.add_vertex(1, x);

  x[0] = 0.0; x[1] = 1.0; x[2] = 0.0;
  editor.add_vertex(2, x);

  x[0] = 0.0; x[1] = 0.0; x[2] = 1.0;
  editor.add_vertex(3, x);

  // Create cells
  editor.init_cells_global(1, 1);
  std::vector<std::size_t> cell_data(4);
  cell_data[0] = 0; cell_data[1] = 1; cell_data[2] = 2; cell_data[3] = 3;
  editor.add_cell(0, cell_data);

  // Close mesh editor
  editor.close();

  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh> MeshFactory::UnitTriangleMesh(MPI_Comm mpi_comm,
                                                    std::string options)
{
  if (MPI::size(mpi_comm) != 1)
    dolfin_error("MeshFactory.cpp",
                 "generate UnitTriangleMesh",
                 "Cannot generate distributed mesh");

  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));

  // Open mesh for editing
  MeshEditor editor;
  editor.open(*mesh, CellType::triangle, 2, 2);

  // Create vertices
  editor.init_vertices_global(3, 3);
  std::vector<double> x(2);
  x[0] = 0.0; x[1] = 0.0;
  editor.add_vertex(0, x);
  x[0] = 1.0; x[1] = 0.0;
  editor.add_vertex(1, x);
  x[0] = 0.0; x[1] = 1.0;
  editor.add_vertex(2, x);

  // Create cells
  editor.init_cells_global(1, 1);
  std::vector<std::size_t> cell_data(3);
  cell_data[0] = 0; cell_data[1] = 1; cell_data[2] = 2;
  editor.add_cell(0, cell_data);

  // Close mesh editor
  editor.close();

  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh>
MeshFactory::UnitHexMesh(MPI_Comm mpi_comm,
                         std::size_t nx, std::size_t ny, std::size_t nz,
                         std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));

  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return mesh;
  }

  MeshEditor editor;
  editor.open(*mesh, CellType::hexahedron, 3, 3);

  // Create vertices and cells:
  editor.init_vertices_global((nx + 1)*(ny + 1)*(nz + 1),
                              (nx + 1)*(ny + 1)*(nz + 1));
  editor.init_cells_global(nx*ny*nz, nx*ny*nz);

  // Storage for vertices
  std::vector<double> x(3);

  const double a = 0.0;
  const double b = 1.0;
  const double c = 0.0;
  const double d = 1.0;
  const double e = 0.0;
  const double f = 1.0;

  // Create main vertices:
  std::size_t vertex = 0;
  for (std::size_t iz = 0; iz <= nz; iz++)
  {
    x[2] = e + ((static_cast<double>(iz))*(f - e)/static_cast<double>(nz));
    for (std::size_t iy = 0; iy <= ny; iy++)
    {
      x[1] = c + ((static_cast<double>(iy))*(d - c)/static_cast<double>(ny));
      for (std::size_t ix = 0; ix <= nx; ix++)
      {
        x[0] = a + ((static_cast<double>(ix))*(b - a)/static_cast<double>(nx));
        editor.add_vertex(vertex, x);
        vertex++;
      }
    }
  }

  // Create cuboids
  std::size_t cell = 0;
  std::vector<std::size_t> v(8);
  for (std::size_t iz = 0; iz < nz; iz++)
    for (std::size_t iy = 0; iy < ny; iy++)
      for (std::size_t ix = 0; ix < nx; ix++)
      {
        v[0] = (iz*(ny + 1) + iy)*(nx + 1) + ix;
        v[1] = v[0] + 1;
        v[2] = v[0] + (nx + 1);
        v[3] = v[1] + (nx + 1);
        v[4] = v[0] + (nx + 1)*(ny + 1);
        v[5] = v[1] + (nx + 1)*(ny + 1);
        v[6] = v[2] + (nx + 1)*(ny + 1);
        v[7] = v[3] + (nx + 1)*(ny + 1);
        editor.add_cell(cell, v);
        ++cell;
      }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh->mpi_comm()))
    MeshPartitioning::build_distributed_mesh(*mesh);

  return mesh;
}
//-----------------------------------------------------------------------------
void MeshFactory::build_rectangle_mesh(std::shared_ptr<Mesh> mesh,
                                       const Point& p0, const Point& p1,
                                       std::size_t nx, std::size_t ny,
                                       std::string diagonal)
{
  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return;
  }

  // Extract minimum and maximum coordinates
  const double x0 = std::min(p0.x(), p1.x());
  const double x1 = std::max(p0.x(), p1.x());
  const double y0 = std::min(p0.y(), p1.y());
  const double y1 = std::max(p0.y(), p1.y());

  const double a = x0;
  const double b = x1;
  const double c = y0;
  const double d = y1;

  if (std::abs(x0 - x1) < DOLFIN_EPS || std::abs(y0 - y1) < DOLFIN_EPS)
  {
    dolfin_error("Rectangle.cpp",
                 "create rectangle",
                 "Rectangle seems to have zero width, height or depth. Consider checking your dimensions");
  }

  if (nx < 1 || ny < 1)
  {
    dolfin_error("RectangleMesh.cpp",
                 "create rectangle",
                 "Rectangle has non-positive number of vertices in some dimension: number of vertices must be at least 1 in each dimension");
  }

  rename("mesh", "Mesh of the unit square (a,b) x (c,d)");
  // Open mesh for editing
  MeshEditor editor;
  editor.open(*mesh, CellType::triangle, 2, 2);

  // Create vertices and cells:
  if (diagonal == "crossed")
  {
    editor.init_vertices_global((nx + 1)*(ny + 1) + nx*ny,
                                  (nx + 1)*(ny + 1) + nx*ny);
    editor.init_cells_global(4*nx*ny, 4*nx*ny);
  }
  else
  {
    editor.init_vertices_global((nx + 1)*(ny + 1), (nx + 1)*(ny + 1));
    editor.init_cells_global(2*nx*ny, 2*nx*ny);
  }

  // Storage for vertices
  std::vector<double> x(2);

  // Create main vertices:
  std::size_t vertex = 0;
  for (std::size_t iy = 0; iy <= ny; iy++)
  {
    x[1] = c + ((static_cast<double>(iy))*(d - c)/static_cast<double>(ny));
    for (std::size_t ix = 0; ix <= nx; ix++)
    {
      x[0] = a + ((static_cast<double>(ix))*(b - a)/static_cast<double>(nx));
      editor.add_vertex(vertex, x);
      vertex++;
    }
  }

  // Create midpoint vertices if the mesh type is crossed
  if (diagonal == "crossed")
  {
    for (std::size_t iy = 0; iy < ny; iy++)
    {
      x[1] = c +(static_cast<double>(iy) + 0.5)*(d - c)/static_cast<double>(ny);
      for (std::size_t ix = 0; ix < nx; ix++)
      {
        x[0] = a + (static_cast<double>(ix) + 0.5)*(b - a)/static_cast<double>(nx);
        editor.add_vertex(vertex, x);
        vertex++;
      }
    }
  }

  // Create triangles
  std::size_t cell = 0;
  if (diagonal == "crossed")
  {
    boost::multi_array<std::size_t, 2> cells(boost::extents[4][3]);
    for (std::size_t iy = 0; iy < ny; iy++)
    {
      for (std::size_t ix = 0; ix < nx; ix++)
      {
        const std::size_t v0 = iy*(nx + 1) + ix;
        const std::size_t v1 = v0 + 1;
        const std::size_t v2 = v0 + (nx + 1);
        const std::size_t v3 = v1 + (nx + 1);
        const std::size_t vmid = (nx + 1)*(ny + 1) + iy*nx + ix;

        // Note that v0 < v1 < v2 < v3 < vmid.
        cells[0][0] = v0; cells[0][1] = v1; cells[0][2] = vmid;
        cells[1][0] = v0; cells[1][1] = v2; cells[1][2] = vmid;
        cells[2][0] = v1; cells[2][1] = v3; cells[2][2] = vmid;
        cells[3][0] = v2; cells[3][1] = v3; cells[3][2] = vmid;

        // Add cells
        for (auto _cell = cells.begin(); _cell != cells.end(); ++_cell)
          editor.add_cell(cell++, *_cell);
      }
    }
  }
  else if (diagonal == "left" || diagonal == "right" || diagonal == "right/left" || diagonal == "left/right")
  {
    std::string local_diagonal = diagonal;
    boost::multi_array<std::size_t, 2> cells(boost::extents[2][3]);
    for (std::size_t iy = 0; iy < ny; iy++)
    {
      // Set up alternating diagonal
      if (diagonal == "right/left")
      {
        if (iy % 2)
          local_diagonal = "right";
        else
          local_diagonal = "left";
      }
      if (diagonal == "left/right")
      {
        if (iy % 2)
          local_diagonal = "left";
        else
          local_diagonal = "right";
      }

      for (std::size_t ix = 0; ix < nx; ix++)
      {
        const std::size_t v0 = iy*(nx + 1) + ix;
        const std::size_t v1 = v0 + 1;
        const std::size_t v2 = v0 + (nx + 1);
        const std::size_t v3 = v1 + (nx + 1);
        std::vector<std::size_t> cell_data;

        if(local_diagonal == "left")
        {
          cells[0][0] = v0; cells[0][1] = v1; cells[0][2] = v2;
          cells[1][0] = v1; cells[1][1] = v2; cells[1][2] = v3;
          if (diagonal == "right/left" || diagonal == "left/right")
            local_diagonal = "right";
        }
        else
        {
          cells[0][0] = v0; cells[0][1] = v1; cells[0][2] = v3;
          cells[1][0] = v0; cells[1][1] = v2; cells[1][2] = v3;
          if (diagonal == "right/left" || diagonal == "left/right")
            local_diagonal = "left";
        }
        editor.add_cell(cell++, cells[0]);
        editor.add_cell(cell++, cells[1]);
      }
    }
  }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return;
  }
}
//-----------------------------------------------------------------------------
void MeshFactory::build_box_mesh(std::shared_ptr<Mesh> mesh,
                                 const Point& p0, const Point& p1,
                                 std::size_t nx, std::size_t ny,
                                 std::size_t nz,
                                 std::string options)
{
  Timer timer("Build BoxMesh");

  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return;
  }

  // Extract minimum and maximum coordinates
  const double x0 = std::min(p0.x(), p1.x());
  const double x1 = std::max(p0.x(), p1.x());
  const double y0 = std::min(p0.y(), p1.y());
  const double y1 = std::max(p0.y(), p1.y());
  const double z0 = std::min(p0.z(), p1.z());
  const double z1 = std::max(p0.z(), p1.z());

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
  editor.open(*mesh, CellType::tetrahedron, 3, 3);

  // Storage for vertex coordinates
  std::vector<double> x(3);

  // Create vertices
  editor.init_vertices_global((nx + 1)*(ny + 1)*(nz + 1),
                              (nx + 1)*(ny + 1)*(nz + 1));
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
  boost::multi_array<std::size_t, 2> cells(boost::extents[6][4]);
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
        for (auto _cell = cells.begin(); _cell != cells.end(); ++_cell)
          editor.add_cell(cell++, *_cell);
      }
    }
  }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh->mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(*mesh);
    return;
  }
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh>
MeshFactory::UnitDiscMesh(MPI_Comm comm, std::size_t n,
                          std::size_t degree, std::size_t gdim,
                          std::string options)
{
  std::shared_ptr<Mesh> mesh(new Mesh(comm));

  if (n == 0)
    dolfin_error("MeshFactory.cpp",
                 "create UnitDiscMesh",
                 "Number of radial divisions must be > 0");
  if (gdim != 2 and gdim != 3)
    dolfin_error("MeshFactory.cpp",
                 "create UnitDiscMesh",
                 "Geometric dimension must be 2 or 3");
  if (degree != 1 and degree != 2)
    dolfin_error("MeshFactory.cpp",
                 "create UnitDiscMesh",
                 "Degree must be 1 or 2");

  MeshEditor editor;
  editor.open(*mesh, 2, gdim, degree);
  editor.init_vertices_global(1 + 3*n*(n + 1),
                              1 + 3*n*(n + 1));

  std::size_t c = 0;
  editor.add_vertex(c, Point(0,0,0));
  ++c;

  for (std::size_t i = 1; i <= n; ++i)
    for (std::size_t j = 0; j < 6*i; ++j)
    {
      double r = (double)i/(double)n;
      double th = 2*M_PI*(double)j/(double)(6*i);
      double x = r*cos(th);
      double y = r*sin(th);
      editor.add_vertex(c, Point(x, y, 0));
      ++c;
    }

  editor.init_cells(6*n*n);

  c = 0;
  std::size_t base_i = 0;
  std::size_t row_i = 1;
  for (std::size_t i = 1; i <= n; ++i)
  {
    std::size_t base_m = base_i;
    base_i = 1 + 3*i*(i - 1);
    std::size_t row_m = row_i;
    row_i = 6*i;

    for (std::size_t k = 0; k != 6; ++k)
      for (std::size_t j = 0; j < (i*2 - 1); ++j)
      {
        std::size_t i0, i1, i2;
        if (j%2 == 0)
        {
          i0 = base_i + (k*i + j/2)%row_i;
          i1 = base_i + (k*i + j/2 + 1)%row_i;
          i2 = base_m + (k*(i-1) + j/2)%row_m;
        }
        else
        {
          i0 = base_m + (k*(i-1) + j/2)%row_m;
          i1 = base_m + (k*(i-1) + j/2 + 1)%row_m;
          i2 = base_i + (k*i + j/2 + 1)%row_i;
        }

        editor.add_cell(c, i0, i1, i2);
        ++c;
      }
  }

  // Initialise entities required for this degree polynomial mesh
  // and allocate space for the point coordinate data

  if (degree == 2)
  {
    editor.init_entities();

    for (EdgeIterator e(*mesh); !e.end(); ++e)
    {
      Point v0 = Vertex(*mesh, e->entities(0)[0]).point();
      Point v1 = Vertex(*mesh, e->entities(0)[1]).point();
      Point pt = e->midpoint();

      if (std::abs(v0.norm() - 1.0) < 1e-6 and
          std::abs(v1.norm() - 1.0) < 1e-6)
        pt *= v0.norm()/pt.norm();

      // Add Edge-based point
      editor.add_entity_point(1, 0, e->index(), pt);
    }
  }

  editor.close();
  return mesh;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Mesh>
MeshFactory::SphericalShellMesh(MPI_Comm mpi_comm, std::size_t degree,
                                std::string options)
{
  if (degree != 1 and degree != 2)
    dolfin_error("MeshFactory.cpp",
                 "generate SphericalShellMesh",
                 "Invalid mesh degree");

  std::shared_ptr<Mesh> mesh(new Mesh(mpi_comm));
  MeshEditor editor;
  const std::size_t tdim = 2;
  const std::size_t gdim = 3;

  editor.open(*mesh, tdim, gdim, degree);

  editor.init_vertices_global(12, 12);

  const double l0 = 2.0/(sqrt(10.0 + 2.0*sqrt(5.0)));
  const double l1 = l0*(1.0 + sqrt(5.0))/2.0;

  // Generate an icosahedron

  editor.add_vertex(0,  Point(  0,  l0, l1));
  editor.add_vertex(1,  Point(  0,  l0, -l1));
  editor.add_vertex(2,  Point(  0, -l0, -l1));
  editor.add_vertex(3,  Point(  0, -l0, l1));
  editor.add_vertex(4,  Point( l1,   0, l0));
  editor.add_vertex(5,  Point(-l1,   0, l0));
  editor.add_vertex(6,  Point(-l1,   0, -l0));
  editor.add_vertex(7,  Point( l1,   0, -l0));
  editor.add_vertex(8,  Point( l0,  l1, 0));
  editor.add_vertex(9,  Point( l0, -l1, 0));
  editor.add_vertex(10, Point(-l0, -l1, 0));
  editor.add_vertex(11, Point(-l0,  l1, 0));

  editor.init_cells_global(20, 20);

  editor.add_cell(0, 0, 4, 8);
  editor.add_cell(1, 0, 5, 11);
  editor.add_cell(2, 1, 6, 11);
  editor.add_cell(3, 1, 7, 8);
  editor.add_cell(4, 2, 6, 10);
  editor.add_cell(5, 2, 7, 9);
  editor.add_cell(6, 3, 4, 9);
  editor.add_cell(7, 3, 5, 10);

  editor.add_cell( 8, 0, 3, 4);
  editor.add_cell( 9, 0, 3, 5);
  editor.add_cell(10, 1, 2, 6);
  editor.add_cell(11, 1, 2, 7);
  editor.add_cell(12, 4, 7, 8);
  editor.add_cell(13, 4, 7, 9);
  editor.add_cell(14, 5, 6, 10);
  editor.add_cell(15, 5, 6, 11);
  editor.add_cell(16, 8, 11, 0);
  editor.add_cell(17, 8, 11, 1);
  editor.add_cell(18, 9, 10, 2);
  editor.add_cell(19, 9, 10, 3);

  if (degree == 2)
  {
    // Initialise entities required for this degree polynomial mesh
    // and allocate space for the point coordinate data
    editor.init_entities();

    for (EdgeIterator e(*mesh); !e.end(); ++e)
    {
      Point v0 = Vertex(*mesh, e->entities(0)[0]).point();
      Point pt = e->midpoint();
      pt *= v0.norm()/pt.norm();

      // Add Edge-based point
      editor.add_entity_point(1, 0, e->index(), pt);
    }
  }

  editor.close();
  return mesh;
}
//-----------------------------------------------------------------------------
