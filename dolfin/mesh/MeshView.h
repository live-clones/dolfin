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

#ifndef __MESHVIEW_H
#define __MESHVIEW_H

#include "Mesh.h"

namespace dolfin
{

  class MeshView : public Mesh
  {
  public:

    // Constructor
    MeshView()
    {};

    /// Constructor from existing Mesh, topological dimension and cell indices
    MeshView(std::shared_ptr<Mesh> mesh,
             std::size_t dim, const std::vector<std::size_t>& indices);

    /// MeshView is a view into a Mesh
    bool is_view() const
    { return true; }

    /// Mesh which this is a view into
    const std::shared_ptr<Mesh> view_mesh() const
    { return mv_mesh; }

    /// Index into original mesh for entities of dimension i
    const std::vector<std::size_t>& view_index(std::size_t i) const
    {
      dolfin_assert(i < mv_index.size());
      return mv_index[i];
    }

    /// Destructor
    ~MeshView()
    {};

  private:

    // Mapping from MeshView indexing to Mesh local index
    // for each topological dimension
    std::vector<std::vector<std::size_t>> mv_index;

    // Mesh which this is a view into
    std::shared_ptr<Mesh> mv_mesh;

  };

}

#endif
