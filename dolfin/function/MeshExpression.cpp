//
// Created by njcs4 on 23/11/17.
//

#include "MeshExpression.h"

using namespace dolfin;

void MeshExpression::eval(Eigen::Ref<Eigen::VectorXd> values,
                  Eigen::Ref<const Eigen::VectorXd> x,
                  const ufc::cell& cell) const
{
  values[0] = (*_mesh_function)[cell.index];
}