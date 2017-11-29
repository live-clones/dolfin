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

#ifndef DOLFIN_MESHEXPRESSION_H
#define DOLFIN_MESHEXPRESSION_H

#include <vector>
#include <ufc.h>
#include <Eigen/Dense>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include "GenericFunction.h"
#include "Expression.h"


namespace dolfin
{
  class MeshExpression : public Expression
  {
  public:
    explicit MeshExpression(std::shared_ptr<MeshFunction<double>> mesh_function) :
        _mesh_function(mesh_function)
    {
      if (mesh_function->dim() < mesh_function->mesh()->topology().dim() - 1)
      {
        dolfin_error("MeshExpression",
                     "Instantiate mesh expression",
                     "MeshFunction topology must be defined on cells or facets");
      }
    }

    explicit MeshExpression(std::shared_ptr<MeshValueCollection<double>> mesh_value_collection) :
        _mesh_value_collection(mesh_value_collection)
    {
      if (mesh_value_collection->dim() < mesh_value_collection->mesh()->topology().dim() - 1)
      {
        dolfin_error("MeshExpression",
                     "Instantiate mesh expression",
                     "MeshValueCollection topology must be defined on cells or facets");
      }
    }

    virtual void eval(Eigen::Ref<Eigen::VectorXd> values,
                      Eigen::Ref<const Eigen::VectorXd> x,
                      const ufc::cell& cell) const;
  private:

    std::shared_ptr<MeshFunction<double>> _mesh_function;
    std::shared_ptr<MeshValueCollection<double>> _mesh_value_collection;
  };
}


#endif
