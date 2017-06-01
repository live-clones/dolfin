// Copyright (C) 2013 Anders Logg
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
// First added:  2013-04-18
// Last changed: 2013-05-30

#ifndef __MESH_POINT_INTERSECTION_H
#define __MESH_POINT_INTERSECTION_H

#include <vector>
#include <memory>

namespace dolfin
{

  // Forward declarations
  class Mesh;
  class Point;

  /// This class represents an intersection between a _Mesh_ and a
  /// _Point_ or a _Mesh_ and a vector defined by two _Point_s.
  /// The resulting intersection is stored as a list of zero
  /// or more _MeshEntity_s.

  class MeshEntityIntersection
  {
  public:

    /// Compute intersection between mesh and point
    MeshEntityIntersection(const Mesh& mesh,
                           const Point& point,
                           const std::size_t t_dim);

    /// Compute intersection between mesh and interval
    MeshEntityIntersection(const Mesh& mesh,
                           const Point& x1,
                           const Point& x2,
                           const std::size_t t_dim);

    ~MeshEntityIntersection() {};

    /// Return the list of (local) indices for intersected entities
    const std::vector<unsigned int>& intersected_entities() const
    { return _intersected_entities; }

  private:

    // The list of (local) indices for intersected entities
    std::vector<unsigned int> _intersected_entities;

  };

}

#endif
