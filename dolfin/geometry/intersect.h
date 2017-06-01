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

#ifndef __INTERSECT_H
#define __INTERSECT_H

#include <memory>

namespace dolfin
{

  // Forward declarations
  class Mesh;
  class Point;
  class MeshPointIntersection;

  /// Compute and return intersection between _Mesh_ and _Point_.
  ///
  /// @param    mesh (_Mesh_)
  ///         The mesh to be intersected.
  /// @param    point (_Point_)
  ///         The point to be intersected.
  ///
  /// @return
  ///     _MeshEntityIntersection_
  ///         The intersection data.
  std::shared_ptr<const MeshEntityIntersection>
  intersect(const Mesh& mesh, const Point& point);

  /// Compute and return intersection between _Mesh_ and _Point_.
  ///
  /// @param    mesh (_Mesh_)
  ///         The mesh to be intersected.
  /// @param    point (_Point_)
  ///         The point to be intersected.
  /// @param    tdim (_size_t_)
  ///         The topology dimension of intersected entities
  ///
  /// @return
  ///     _MeshEntityIntersection_
  ///         The intersection data.
  std::shared_ptr<const MeshEntityIntersection>
  intersect(const Mesh& mesh, const Point& point,
            const std::size_t tdim);

  /// Compute and return intersection between _Mesh_ and vector between
  /// two _Point_s.
  ///
  /// @param    mesh (_Mesh_)
  ///         The mesh to be intersected.
  /// @param    x1 (_Point_)
  ///         The vector origin.
  /// @param    x2 (_Point_)
  ///         The vector destination.
  ///
  /// @return
  ///     _MeshEntityIntersection_
  ///         The intersection data.
  std::shared_ptr<const MeshEntityIntersection>
  intersect(const Mesh& mesh, const Point& x1, const Point& x2);


  /// Compute and return intersection between _Mesh_ and vector between
  /// two _Point_s.
  ///
  /// @param    mesh (_Mesh_)
  ///         The mesh to be intersected.
  /// @param    x1 (_Point_)
  ///         The vector origin.
  /// @param    x2 (_Point_)
  ///         The vector destination.
  /// @param    tdim (_size_t_)
  ///         The topology dimension of intersected entities
  ///
  /// @return
  ///     _MeshEntityIntersection_
  ///         The intersection data.
  std::shared_ptr<const MeshEntityIntersection>
  intersect(const Mesh& mesh, const Point& x1, const Point& x2,
            const std::size_t tdim);

}

#endif
