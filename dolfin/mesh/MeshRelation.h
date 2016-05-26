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
//

#ifndef __MESH_RELATION_H
#define __MESH_RELATION_H

#include <vector>
#include <memory>

#include "MeshHierarchy.h"

namespace dolfin
{
  class Mesh;

  /// MeshRelation encapsulates the relationships which may exist between two Meshes
  /// or which may exist in a Mesh as a result of being related to another Mesh

  class MeshRelation
  {
  public:
    /// Constructor
    MeshRelation()
    {}

    /// Constructor with source and destination meshes
    MeshRelation(std::shared_ptr<const Mesh> m1,
                 std::shared_ptr<const Mesh> m2)
      : src_mesh(m1), dest_mesh(m2)
    {}

    /// Destructor
    ~MeshRelation()
    {}

    /// Get the source mesh
    std::shared_ptr<const Mesh> source_mesh() const
    { return src_mesh; }

    /// Get the destination mesh
    std::shared_ptr<const Mesh> destination_mesh() const
    { return dest_mesh; }

    /// Initialise the mapping of entities of dim from the source to destination mesh
    void init(unsigned int dim, std::vector<std::size_t> map_back_to_src);

    /// Initialise the mapping of entities of dim from the source to destination mesh
    /// when there are multiple mappings (e.g. parent -> child cells)
    void init(unsigned int dim, std::vector<std::vector<std::size_t>> map_back_to_src);

    /// Get the number of entities of dimension dim, pointed to by entity index in destination mesh
    /// e.g. if destination mesh is a coarse mesh, and source mesh is fine, each coarse cell
    /// could point to multiple fine cells.
    unsigned int num_entities(unsigned int dim, unsigned int index);

    /// Pointer to the indices on the source mesh, for index in the destination mesh.
    /// In many cases, there will just be one index. For example, vertices are always mapped one-to-one, so this
    /// will just be a map back from the index in the destination mesh to the index of the same vertex in the
    /// source mesh.
    unsigned int* entities(unsigned int dim, unsigned int index);

    // Get mapping to single entity on source mesh from destination mesh
    unsigned int entity(unsigned int dim, unsigned int index)
    {
      dolfin_assert(dest_to_src.size() > dim);
      dolfin_assert(dest_to_src_offset[dim].empty());
      auto& entity_indices = dest_to_src[dim];
      dolfin_assert(index < entity_indices.size());
      return entity_indices[index];
    }

    /// Get size of relationship map back to source mesh. Should be same size as number of entities of
    /// dimension dim in destination mesh, if calculated.
    unsigned int size(unsigned int dim)
    {
      dolfin_assert(dim < dest_to_src.size());
      return dest_to_src[dim].size();
    }

  private:

    // FIXME: remove all this stuff about edge_to_global_vertex from here
    friend class MeshHierarchy;
    friend class PlazaRefinementND;
    // Map from edge of parent Mesh to new vertex in child Mesh
    // as calculated during ParallelRefinement process
    std::shared_ptr<const std::map<std::size_t, std::size_t> > edge_to_global_vertex;

    // Shared pointers to source and destination meshes
    std::shared_ptr<const Mesh> src_mesh;
    std::shared_ptr<const Mesh> dest_mesh;

    // Mapping of entities of each dimension from the
    // source to destination mesh. In the case of
    // a coarse->fine mesh relationship, each coarse
    // cell may map to multiple fine cells.
    // Mappings of all entity dims may not be computed.
    std::vector<std::vector<unsigned int>> dest_to_src;

    // Offset vectors, only used when mapping to multiple entities
    // on the destination mesh, e.g. parent->children
    std::vector<std::vector<unsigned int>> dest_to_src_offset;

  };
}

#endif
