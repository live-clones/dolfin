// Copyright (C) 2016 Chris Richardson, Garth N. Wells and Cian Wilson
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
// Modified by Cian Wilson, 2016

#ifndef __COORDINATE_DOF_MAP_H
#define __COORDINATE_DOF_MAP_H

#include <cstdlib>
#include <vector>

namespace dolfin
{

  class Mesh;

  /// This class handles the mapping of coordinate degrees of freedom.

  class CoordinateDofMap
  {
  public:

    /// Create coordinate dof map
    CoordinateDofMap();

    /// Destructor
    ~CoordinateDofMap();

    void init(const Mesh& mesh);

    /// Return true if the dofmap is empty
    bool empty() const
    { return _dofmap.empty(); }

    /// Return array of connections for given cell
    const unsigned int* operator() (std::size_t cell) const
    {
      return &_dofmap[cell*_dofs_per_cell];
    }

    /// Return the number of dofs reuiqred for cell
    std::size_t num_cell_dofs() const
    { return _dofs_per_cell; }

  private:

    // Cell-local-to-dof map (dofs for cell dofmap[i])
    std::vector<unsigned int> _dofmap;

    // Number of dofs specified per cell (fixed for now)
    std::size_t _dofs_per_cell;

  };
}

#endif
