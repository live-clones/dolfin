// Copyright (C) 2016 Chris Richardson
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

#ifndef __DOLFIN_MESHSLICER_H
#define __DOLFIN_MESHSLICER_H

#include <boost/multi_array.hpp>

namespace dolfin
{

  class Plane;

  /// Implements the sub-division of a Cell or Mesh by a surface

  class MeshSlicer
  {
  public:

    /// Return a new Mesh, adding all points to the given Mesh
    /// which are on the Plane where it intersects the Mesh
    static Mesh add(const Mesh& mesh,
                    const Plane& plane);

    /// Return a 2D Mesh of triangles representing the cut
    /// of the Mesh with the Plane
    static Mesh cut(const Mesh& mesh,
                    const Plane& plane);

  private:

    // Generate volume or surface mesh
    static void create_mesh(Mesh& mesh,
                            const boost::multi_array<double, 2>& geometry,
                            const std::vector<std::size_t>& topology,
                            std::size_t tdim);

  };

}

#endif
