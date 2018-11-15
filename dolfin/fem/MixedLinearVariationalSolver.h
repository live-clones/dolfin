// Copyright (C) 2008-2011 Anders Logg and Garth N. Wells
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
// Modified by Cecile Daversin-Catty, 2017.
//
// First added:  2017-07-21
// Last changed: 2017-07-21

#ifndef __MIXED_LINEAR_VARIATIONAL_SOLVER_H
#define __MIXED_LINEAR_VARIATIONAL_SOLVER_H

#include <dolfin/common/Variable.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/la/KrylovSolver.h>

namespace dolfin
{

  // Forward declarations
  class MixedLinearVariationalProblem;
  class PETScNestMatrix;

  /// This class implements a solver for mixed linear variational problems.

  class MixedLinearVariationalSolver : public Variable
  {
  public:
    typedef std::tuple<std::vector<std::shared_ptr<GenericMatrix>>,
      std::vector<std::shared_ptr<GenericVector>>,
      std::vector<std::shared_ptr<GenericVector>> > assembled_system_type;

    /// Create linear variational solver for given problem
    explicit MixedLinearVariationalSolver();
    explicit MixedLinearVariationalSolver(std::shared_ptr<MixedLinearVariationalProblem> problem);

    /// Block-by-block assembly
    assembled_system_type assemble_system();

    /// Solve variational problem
    void solve();
    void solve(assembled_system_type assembled_system);

    /// Default parameter values
    static Parameters default_parameters()
    {
      Parameters p("linear_variational_solver");

      p.add("linear_solver", "default");
      p.add("preconditioner", "default");
      p.add("symmetric", false);

      p.add("print_rhs", false);
      p.add("print_matrix", false);
      p.add(LUSolver::default_parameters());
      p.add(KrylovSolver::default_parameters());

      return p;
    }

  private:

    // The linear problem
    std::shared_ptr<MixedLinearVariationalProblem> _problem;

  };

}

#endif
