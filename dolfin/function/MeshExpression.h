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

  /// The MeshExpression class is an Expression which yields piecewise
  /// constant values defined on the entity topology of a MeshFunction
  /// or MeshValueCollection.
  /// This is useful, for example, in the case that a boundary condition
  /// has facet-dependent values defined on each facet in the mesh. Or,
  /// for example, in the case that each cell in a mesh has piecewise
  /// constant material data associated with it. These values may be
  /// straightforwardly read in from file via the MeshFunction or
  /// MeshValueCollection interface.

  class MeshExpression : public Expression
  {
  public:

    /// Create a MeshExpression where the facet or cell values are stored in
    /// a (_MeshFunction_)
    ///
    /// @param    mesh_function (_MeshFunction_)
    ///         The entity values.
    explicit MeshExpression(std::shared_ptr<MeshFunction<double>> mesh_function);

    /// Create a MeshExpression where the facet or cell values are stored in
    /// a (_MeshValueCollection_)
    ///
    /// @param    mesh_value_collection (_MeshValueCollection_)
    ///         The entity values.
    explicit MeshExpression(std::shared_ptr<MeshValueCollection<double>> mesh_value_collection);

    virtual void eval(Eigen::Ref<Eigen::VectorXd> values,
                      Eigen::Ref<const Eigen::VectorXd> x,
                      const ufc::cell& cell) const
    {
      eval_function(values, x, cell);
    }

  private:

    std::function<void(Eigen::Ref<Eigen::VectorXd> values,
                       Eigen::Ref<const Eigen::VectorXd> x,
                       const ufc::cell& cell)> eval_function;

    void eval_mesh_function(Eigen::Ref<Eigen::VectorXd> values,
                            Eigen::Ref<const Eigen::VectorXd> x,
                            const ufc::cell& cell) const;

    void eval_mesh_value_collection(Eigen::Ref<Eigen::VectorXd> values,
                                    Eigen::Ref<const Eigen::VectorXd> x,
                                    const ufc::cell& cell) const;

    std::shared_ptr<MeshFunction<double>> _mesh_function;
    std::shared_ptr<MeshValueCollection<double>> _mesh_value_collection;
  };
}


#endif
