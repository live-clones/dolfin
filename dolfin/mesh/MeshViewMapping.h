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

#ifndef __MESH_VIEW_MAPPING_H
#define __MESH_VIEW_MAPPING_H

namespace dolfin
{

  class MeshViewMapping
  {
  public:

    /// Default Constructor
    MeshViewMapping() : _mesh(0)
    {}

    /// Create a new Mesh based on the Meshfunction marker, where it has a value equal
    /// to tag.
    static Mesh create_from_marker(const MeshFunction<std::size_t>& marker, std::size_t tag);

  private:

    // The mesh which this mapping points to
    std::shared_ptr<const Mesh> _mesh;

    // Map to cells in _mesh
    std::vector<std::size_t> cell_map;

  };

}

#endif
