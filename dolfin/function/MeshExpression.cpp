// Copyright (C) 2017 Nate Sime
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

#include <dolfin/mesh/Cell.h>

#include "MeshExpression.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshExpression::eval(Eigen::Ref<Eigen::VectorXd> values,
                  Eigen::Ref<const Eigen::VectorXd> x,
                  const ufc::cell& cell) const
{
  dolfin_assert(_mesh_function);
  dolfin_assert(_mesh_function->mesh());
  const std::size_t mesh_tdim = _mesh_function->mesh()->topology().dim();

  dolfin_assert(_mesh_function->dim() > mesh_tdim - 1);

  if ((cell.local_facet < 0) and ((mesh_tdim - 1) == _mesh_function->dim()))
  {
    dolfin_error("MeshExpression::eval",
                 "evaluate facet MeshExpression on a cell",
                 "MeshExpression topological dimension permits solely facet integration.");
  }

  if ((cell.local_facet < 0) or (mesh_tdim == _mesh_function->dim()))
  {
    values[0] = (*_mesh_function)[cell.index];
  }
  else
  {
    const unsigned int facet_idx =
        Cell(*_mesh_function->mesh(), cell.index).entities(_mesh_function->dim())[cell.local_facet];
    values[0] = (*_mesh_function)[facet_idx];
  }
}
//-----------------------------------------------------------------------------