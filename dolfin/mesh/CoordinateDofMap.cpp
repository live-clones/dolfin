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

#include "CoordinateDofMap.h"
#include "Mesh.h"
#include "CellType.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
CoordinateDofMap::CoordinateDofMap()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CoordinateDofMap::~CoordinateDofMap()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void CoordinateDofMap::init(const Mesh& mesh)
{
  std::size_t tdim = mesh.topology().dim();

  // Initialise with vertices only (for now)
  _dofs_per_cell = mesh.type().num_vertices();

  const MeshConnectivity& connectivity = mesh.topology()(tdim, 0);
  _dofmap = connectivity();
}
//-----------------------------------------------------------------------------
