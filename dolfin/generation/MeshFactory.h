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

#ifndef __MESHFACTORY_H
#define __MESHFACTORY_H

#include <memory>

#include <dolfin/common/MPI.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

  class MeshFactory
  {
  public:

    /// Create a uniform finite element _Mesh_ over the unit cube
    /// [0,1] x [0,1] x [0,1].
    ///
    /// *Arguments*
    ///     comm (MPI_Comm)
    ///         MPI communicator
    ///     nx (std::size_t)
    ///         Number of cells in :math:`x` direction.
    ///     ny (std::size_t)
    ///         Number of cells in :math:`y` direction.
    ///     nz (std::size_t)
    ///         Number of cells in :math:`z` direction.
    ///     options (std::string)
    ///         Any options for mesh creation
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         std::shared_ptr<Mesh> mesh = MeshFactorty::UnitCubeMesh(MPI_COMM_WORLD, 32, 32, 32);
    ///
    static std::shared_ptr<Mesh>
      UnitCubeMesh(MPI_Comm comm, std::size_t nx, std::size_t ny,
                   std::size_t nz, std::string options="");

    /// Create a uniform finite element _Mesh_ over the rectangular
    /// prism spanned by the two _Point_s p0 and p1. The order of the
    /// two points is not important in terms of minimum and maximum
    /// coordinates.
    ///
    /// *Arguments*
    ///     comm (MPI_Comm)
    ///         MPI communicator
    ///     p0 (_Point_)
    ///         First point.
    ///     p1 (_Point_)
    ///         Second point.
    ///     nx (double)
    ///         Number of cells in :math:`x`-direction.
    ///     ny (double)
    ///         Number of cells in :math:`y`-direction.
    ///     nz (double)
    ///         Number of cells in :math:`z`-direction.
    ///     options (std::string)
    ///         Any options for mesh creation
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         // Mesh with 8 cells in each direction on the
    ///         // set [-1,2] x [-1,2] x [-1,2].
    ///         Point p0(-1, -1, -1);
    ///         Point p1(2, 2, 2);
    ///         mesh = MeshFactory::BoxMesh(MPI_COMM_WORLD, p0, p1, 8, 8, 8);
    ///
    static std::shared_ptr<Mesh>
      BoxMesh(MPI_Comm comm,
              const Point& p0, const Point& p1,
              std::size_t nx, std::size_t ny, std::size_t nz,
              std::string options="");

    /// Create a uniform finite element _Mesh_ over the unit square
    /// [0,1] x [0,1].
    ///
    /// *Arguments*
    ///     mpi_comm (MPI_Comm)
    ///         MPI communicator
    ///     nx (std::size_t)
    ///         Number of cells in horizontal direction.
    ///     ny (std::size_t)
    ///         Number of cells in vertical direction.
    ///     options (std::string)
    ///         Options argument: 'left', 'right', or 'crossed' indicating
    ///         the direction of the diagonals.
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         std::shared_ptr<Mesh> mesh1 = MeshFactory::UnitSquareMesh(MPI_COMM_WORLD, 32, 32);
    ///         std::shared_ptr<Mesh> mesh2 = MeshFactory::UnitSquareMesh(MPI_COMM_WORLD, 32, 32, "crossed");
    ///
    static std::shared_ptr<Mesh>
      UnitSquareMesh(MPI_Comm mpi_comm, std::size_t nx, std::size_t ny,
                     std::string options="right");

    /// *Arguments*
    ///     comm (MPI_Comm)
    ///         MPI communicator
    ///     p0 (_Point_)
    ///         First point.
    ///     p1 (_Point_)
    ///         Second point.
    ///     nx (double)
    ///         Number of cells in :math:`x`-direction.
    ///     ny (double)
    ///         Number of cells in :math:`y`-direction.
    ///     options (std::string)
    ///         Options argument: 'left', 'right', or 'crossed' indicating
    ///         the direction of the diagonals.
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         // Mesh with 8 cells in each direction on the
    ///         // set [-1,2] x [-1,2]
    ///         Point p0(-1, -1);
    ///         Point p1(2, 2);
    ///         std::shared_ptr<Mesh> mesh = MeshFactory::RectangleMesh(MPI_COMM_WORLD, p0, p1, 8, 8);
    ///
    static std::shared_ptr<Mesh>
      RectangleMesh(MPI_Comm mpi_comm, const Point& p0, const Point& p1,
                    std::size_t nx, std::size_t ny,
                    std::string options="right");

    /// Create a uniform finite element _Mesh_ over the unit interval
    /// [0,1].
    ///
    /// *Arguments*
    ///     mpi_comm (MPI_Comm)
    ///         MPI communicator
    ///     nx (std::size_t)
    ///         Number of cells in horizontal direction.
    ///     options (std::string)
    ///         Options argument.
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         std::shared_ptr<Mesh> mesh1 = MeshFactory::UnitIntervalMesh(MPI_COMM_WORLD, 32);
    ///
    static std::shared_ptr<Mesh>
      UnitIntervalMesh(MPI_Comm mpi_comm, std::size_t nx,
                       std::string options="")
    {
      return MeshFactory::IntervalMesh(mpi_comm, nx, 0.0, 1.0, options);
    }

    /// *Arguments*
    ///     comm (MPI_Comm)
    ///         MPI communicator
    ///     p0 (_Point_)
    ///         First point.
    ///     p1 (_Point_)
    ///         Second point.
    ///     nx (double)
    ///         Number of cells in :math:`x`-direction.
    ///     ny (double)
    ///         Number of cells in :math:`y`-direction.
    ///     options (std::string)
    ///         Options argument.
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         double a = 1, b = 2;
    ///         std::shared_ptr<Mesh> mesh = MeshFactory::IntervalMesh(MPI_COMM_WORLD, a, b, 8);
    ///
    static std::shared_ptr<Mesh>
      IntervalMesh(MPI_Comm mpi_comm, std::size_t nx, double a, double b,
                   std::string options="");

    /// A mesh consisting of a single tetrahedron with vertices at
    ///
    ///   (0, 0, 0)
    ///   (1, 0, 0)
    ///   (0, 1, 0)
    ///   (0, 0, 1)
    /// Useful for testing
    static std::shared_ptr<Mesh>
      UnitTetrahedronMesh(MPI_Comm mpi_comm, std::string options="");


    /// A mesh consisting of a single triangle with vertices at
    ///
    ///   (0, 0)
    ///   (1, 0)
    ///   (0, 1)
    /// Useful for testing
    static std::shared_ptr<Mesh>
      UnitTriangleMesh(MPI_Comm mpi_comm, std::string options="");

    /// EXPERIMENTAL
    /// A hexahedral cubic box mesh, for testing purposes, not usable
    /// with anything else
    static std::shared_ptr<Mesh>
      UnitHexMesh(MPI_Comm mpi_comm, std::size_t nx, std::size_t ny,
                  std::size_t nz, std::string options="");

    /// EXPERIMENTAL
    /// A quadrilateral square mesh, for testing purposes, not usable
    /// with anything else
    static std::shared_ptr<Mesh>
      UnitQuadMesh(MPI_Comm mpi_comm, std::size_t nx, std::size_t ny,
                   std::string options="");

    /// EXPERIMENTAL
    /// A mesh of the Unit Circle of radius 1, centered at (0,0,0)
    /// either in a 2D or 3D geometry, with a linear or quadratic mesh degree.
    static std::shared_ptr<Mesh>
      UnitDiscMesh(MPI_Comm mpi_comm, std::size_t n,
                   std::size_t degree, std::size_t gdim,
                   std::string options="");

    /// EXPERIMENTAL
    /// An icosahedral mesh, approximating a spherical shell, optionally in
    /// linear or quadratic degree. Can be refined to get a finer mesh.
    static std::shared_ptr<Mesh>
      SphericalShellMesh(MPI_Comm mpi_comm, std::size_t degree,
                         std::string options="");

  private:
    // Generate a rectangle mesh of size nx*ny between points p0 and p1
    static void build_rectangle_mesh(std::shared_ptr<Mesh> mesh,
                                     const Point& p0, const Point& p1,
                                     std::size_t nx, std::size_t ny,
                                     std::string options);

    // Generate a box mesh of size nx*ny*nz between points p0 and p1
    static void build_box_mesh(std::shared_ptr<Mesh> mesh,
                               const Point& p0, const Point& p1,
                               std::size_t nx, std::size_t ny, std::size_t nz,
                               std::string options);

  };
}

#endif
