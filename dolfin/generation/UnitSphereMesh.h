// Copyright (C) 2017 Chris Richardson
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

#ifndef __UNIT_SPHERE_MESH_H
#define __UNIT_SPHERE_MESH_H

#include <dolfin/common/MPI.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{
  /// Icosahedral mesh, with degree=1
  /// or degree=2 with a refinement option.

  class UnitSphereMesh
  {
  public:

    /// Create a spherical mesh for testing
    static Mesh create(MPI_Comm comm, std::size_t nrefine, std::size_t degree)
    {
      Mesh mesh(comm);
      build(mesh, nrefine, degree);
      return mesh;
    }

  private:

    static void build(Mesh& mesh, std::size_t nrefine, std::size_t degree);

    // Move points to a spherical surface
    static void move_surface_points(Mesh& mesh);

  };

}

#endif
