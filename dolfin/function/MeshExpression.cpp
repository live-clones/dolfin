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

#include "MeshExpression.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshExpression::eval(Eigen::Ref<Eigen::VectorXd> values,
                  Eigen::Ref<const Eigen::VectorXd> x,
                  const ufc::cell& cell) const
{
  info("ufc cell dim %d, index %d, local facet %d", cell.topological_dimension, cell.index, cell.local_facet);
//  dolfin_assert(cell.topological_dimension == _mesh_function->dim());
  values[0] = (*_mesh_function)[cell.index];
}
//-----------------------------------------------------------------------------