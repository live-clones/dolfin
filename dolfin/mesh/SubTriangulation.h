// Copyright (C) 2013 Melior Technology
//
// This file is part of Narwal.
//
// Narwal is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Narwal is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Narwal. If not, see <http://www.gnu.org/licenses/>.

#ifndef __DOLFIN_SUBTRIANGULATION_H
#define __DOLFIN_SUBTRIANGULATION_H

#include <boost/multi_array.hpp>

namespace dolfin
{

  class Plane;

  /// Implements the sub-division of a Cell or Mesh by a surface

  class SubTriangulation
  {
  public:

    /// Build a SubTriangulation on a Cell, as it is intersected
    /// by the given plane
    SubTriangulation(const Cell& cell,
                     const Plane& plane);

    /// Build a SubTriangulation on an entire Mesh, based on
    /// the intersection with the given plane
    SubTriangulation(const Mesh& mesh,
                     const Plane& plane);

    /// Destructor
    ~SubTriangulation() {}

    /// Small Mesh which is formed by cutting the divided cell into
    /// subcells
    const Mesh& mesh() const
    { return _mesh; }

  private:

    // Volume mesh including all intersection points with surface
    Mesh _mesh;

    // Rebuild subtriangulation with new surface
    void build(const Plane& surface);

    // Generate volume or surface mesh
    void create_mesh(Mesh& mesh,
                     const boost::multi_array<double, 2>& geometry,
                     const std::vector<std::size_t>& topology,
                     std::size_t tdim);

  };

}


#endif
