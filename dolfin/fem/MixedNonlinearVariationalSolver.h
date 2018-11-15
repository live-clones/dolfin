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
// Modified by Cecile Daversin-Catty, 2018.
//
// First added:  2017-10-03
// Last changed: 2017-10-03

#ifndef __MIXED_NONLINEAR_VARIATIONAL_SOLVER_H
#define __MIXED_NONLINEAR_VARIATIONAL_SOLVER_H

#include <dolfin/nls/NonlinearProblem.h>
#include <dolfin/nls/NewtonSolver.h>
#include <dolfin/nls/PETScSNESSolver.h>
#include "MixedNonlinearVariationalProblem.h"
#include "SystemAssembler.h"

namespace dolfin
{

  // Forward declarations
  class MixedNonlinearVariationalProblem;
  class PETScNestMatrix;

  /// This class implements a solver for mixed nonlinear variational problems.

  class MixedNonlinearVariationalSolver : public Variable
  {
  public:
    /* typedef std::tuple<std::vector<std::shared_ptr<GenericVector>>, */
    /*   std::vector<std::shared_ptr<GenericMatrix>> > assembled_system_type; */

    /// Create linear variational solver for given problem
    explicit MixedNonlinearVariationalSolver();
    explicit MixedNonlinearVariationalSolver(std::shared_ptr<MixedNonlinearVariationalProblem> problem);
    /// Solve variational problem
    std::pair<std::size_t, bool> solve();

    /// Default parameter values
    static Parameters default_parameters()
    {
      Parameters p("nonlinear_variational_solver");

      p.add("symmetric", false);
      p.add("print_rhs", false);
      p.add("print_matrix", false);

      std::set<std::string> nonlinear_solvers = {"newton"};
      std::string default_nonlinear_solver = "newton";
      p.add(NewtonSolver::default_parameters());

      #ifdef HAS_PETSC
      p.add(PETScSNESSolver::default_parameters());
      nonlinear_solvers.insert("snes");
      #endif

      p.add("nonlinear_solver", default_nonlinear_solver, nonlinear_solvers);

      return p;
    }

  private:

    // Nonlinear (algebraic) problem
    class MixedNonlinearDiscreteProblem : public NonlinearProblem
    {
    public:

      // Constructor
      MixedNonlinearDiscreteProblem(
        std::shared_ptr<const MixedNonlinearVariationalProblem> problem,
        std::shared_ptr<const MixedNonlinearVariationalSolver> solver);

      // Destructor
      ~MixedNonlinearDiscreteProblem();

      // Compute F at current point x
      virtual void F(GenericVector& b, const GenericVector& x);

      // Compute J = F' at current point x
      virtual void J(GenericMatrix& A, const GenericVector& x);

      std::vector<std::shared_ptr<GenericMatrix>> _Js;

    private:

      // Problem and solver objects
      std::shared_ptr<const MixedNonlinearVariationalProblem> _problem;
      std::shared_ptr<const MixedNonlinearVariationalSolver> _solver;

    };

    // The nonlinear problem
    std::shared_ptr<MixedNonlinearVariationalProblem> _problem;

    // The nonlinear discrete problem
    std::shared_ptr<MixedNonlinearDiscreteProblem> nonlinear_problem;

    // The Newton solver
    std::shared_ptr<NewtonSolver> newton_solver;

    #ifdef HAS_PETSC
    // Or, alternatively, the SNES solver
    std::shared_ptr<PETScSNESSolver> snes_solver;
    #endif

  };
}

#endif
