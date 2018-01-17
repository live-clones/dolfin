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

#ifndef __MESH_VIEW_H
#define __MESH_VIEW_H

namespace dolfin
{
  // Forward declarations
  class Mesh;
  template <typename T> class MeshFunction;

  class MeshView
  {
    /// Mapping from a child Mesh to a parent Mesh for vertices
    /// and cells, using local indexing

  public:

    /// Constructor
    MeshView(std::shared_ptr<const Mesh> parent_mesh,
             std::vector<std::size_t>& vertex_map,
             std::vector<std::size_t>& cell_map)
      : _mesh(parent_mesh), _vertex_map(vertex_map), _cell_map(cell_map)
    {
      // Do nothing
    }

    /// Access parent mesh
    std::shared_ptr<const Mesh> mesh() const
    {
      return _mesh;
    }

    /// Map to vertices of parent mesh
    const std::vector<std::size_t>& vertex_map() const
    {
      return _vertex_map;
    }

    /// Map to cells of parent mesh
    const std::vector<std::size_t>& cell_map() const
    {
      return _cell_map;
    }

    /// Create a new Mesh based on the Meshfunction marker, where it has a value equal
    /// to tag, setting the MeshViewMapping in MeshTopology accordingly.
    /// FIXME: this could be a free function
    static Mesh create(const MeshFunction<std::size_t>& marker, std::size_t tag);

  private:

    // The parent mesh which this mapping points to
    std::shared_ptr<const Mesh> _mesh;

    // Map to vertices in _mesh
    std::vector<std::size_t> _vertex_map;

    // Map to cells in _mesh
    std::vector<std::size_t> _cell_map;

  };

}

#endif
