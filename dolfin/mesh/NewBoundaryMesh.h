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
//

#ifndef __NEW_BOUNDARY_MESH_H
#define __NEW_BOUNDARY_MESH_H

#include "Mesh.h"
#include "MeshRelation.h"

namespace dolfin
{
  class NewBoundaryMesh
  {
  public:

    /// Return a Mesh of the external boundary
    static std::shared_ptr<const Mesh> create_boundary_mesh(std::shared_ptr<const Mesh> mesh)
    { return create_boundary_relation(mesh)->destination_mesh(); }

    /// Return a MeshRelation of the Mesh with its external boundary (the boundary mesh is available as
    /// MeshRelation.destination_mesh()).
    static std::shared_ptr<const MeshRelation> create_boundary_relation(std::shared_ptr<const Mesh> mesh);

    /// Return a MeshRelation of the Mesh with a new Mesh described by the entity markers on dimension tdim.
    /// The new mesh will be available as MeshRelation.destination_mesh().
    static std::shared_ptr<const MeshRelation> create_relation(std::shared_ptr<const Mesh> mesh,
							       std::vector<std::size_t> markers,
							       std::size_t tdim);

    /// Initialise the relationship between edges on the new Mesh and edges on the original mesh
    /// FIXME: maybe this should be a private method, called by option from create_relation
    static void initialise_edge_relation(std::shared_ptr<MeshRelation> meshrelation);


  };

}

#endif
