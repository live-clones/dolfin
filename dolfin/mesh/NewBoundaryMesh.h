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

namespace dolfin
{
  class MeshRelation;

  class NewBoundaryMesh
  {
  public:

    static std::shared_ptr<const MeshRelation> create_boundary(std::shared_ptr<const Mesh> mesh);

    static std::shared_ptr<const MeshRelation> create(std::shared_ptr<const Mesh> mesh,
                                                      std::vector<std::size_t> markers,
                                                      std::size_t tdim);

  };

}

#endif
