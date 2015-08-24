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
  enum class MeshOptions : int { none = 0, left = 1, right = 2, crossed = 3, alternating = 4 };
  inline MeshOptions operator&(MeshOptions a, MeshOptions b)
  { return static_cast<MeshOptions>(static_cast<int>(a)&static_cast<int>(b));}
  inline MeshOptions operator|(MeshOptions a, MeshOptions b)
  { return static_cast<MeshOptions>(static_cast<int>(a)|static_cast<int>(b));}

  class MeshFactory
  {
  public:
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
    ///     options
    ///         Options argument: 'left', 'right', or 'crossed' indicating
    ///         the direction of the diagonals.
    ///
    /// *Example*
    ///     .. code-block:: c++
    ///
    ///         UnitSquareMesh mesh1(MPI_COMM_WORLD, 32, 32);
    ///         UnitSquareMesh mesh2(MPI_COMM_WORLD, 32, 32, MeshOptions::crossed);
    ///

    static std::shared_ptr<Mesh>
      UnitSquareMesh(MPI_Comm mpi_comm, std::size_t nx, std::size_t ny,
                               MeshOptions options=MeshOptions::left);

    static std::shared_ptr<Mesh>
      RectangleMesh(MPI_Comm mpi_comm, Point p0, Point p1,
                    std::size_t nx, std::size_t ny,
                    MeshOptions options=MeshOptions::left);

  private:
    // Generate a rectangle mesh of size nx*ny between points p0 and p1
    static void build_rectangle_mesh(std::shared_ptr<Mesh> mesh,
                                     Point p0, Point p1,
                                     std::size_t nx, std::size_t ny,
                                     MeshOptions options);

  };
}

#endif
