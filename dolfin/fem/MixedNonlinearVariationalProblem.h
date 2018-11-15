// Copyright (C) 2011 Anders Logg
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
// First added:  2018-10-03
// Last changed: 2018-10-03

#ifndef __MIXED_NONLINEAR_VARIATIONAL_PROBLEM_H
#define __MIXED_NONLINEAR_VARIATIONAL_PROBLEM_H

#include <memory>
#include <vector>
#include <dolfin/common/Hierarchical.h>

namespace dolfin
{

  // Forward declarations
  class DirichletBC;
  class Form;
  class Function;
  class FunctionSpace;

  /// This class represents a mixed nonlinear variational problem:
  ///
  /// Find u = (u_1, ..., u_n) in V = V1 x ... x Vn such that
  ///
  ///     F(u; v) = 0  for all v in V^,
  ///
  /// where V is the trial space and V^ is the test space.

  class MixedNonlinearVariationalProblem : public Hierarchical<MixedNonlinearVariationalProblem>
  {
  public:
    /// Create mixed nonlinear variational problem
    /// with a list of boundary conditions
    typedef std::vector<std::vector<std::shared_ptr<const Form>>> form_list_type;
    MixedNonlinearVariationalProblem(form_list_type F,
				     std::vector<std::shared_ptr<Function>> u,
				     std::vector<std::shared_ptr<const DirichletBC>> bcs,
				     form_list_type J={{nullptr}});

    /// Set the bounds for bound constrained solver
    void set_bounds(const std::vector<Function>& lb_func,
		    const std::vector<Function>& ub_func);

    /// Set the bounds for bound constrained solver
    void set_bounds(std::vector<std::shared_ptr<const GenericVector>> lb,
                    std::vector<std::shared_ptr<const GenericVector>> ub);

    /// Return residual form
    form_list_type residual_form() const;
    std::shared_ptr<const Form> residual_form(int i, int j=0) const;
    
    /// Return Jacobian form
    form_list_type jacobian_form() const;
    std::shared_ptr<const Form> jacobian_form(int i, int j=0) const;

    /// Return solution variable
    std::vector<std::shared_ptr<Function>> solution();
    std::shared_ptr<Function> solution(int i);
    /// Return solution variable - const version
    //std::vector<std::shared_ptr<const Function>> solution() const;
    const std::vector<std::shared_ptr<Function>> solution() const;
    std::shared_ptr<const Function> solution(int i) const;

    /// Return boundary conditions
    std::vector<std::vector<std::shared_ptr<const DirichletBC>>> bcs() const;
    std::vector<std::shared_ptr<const DirichletBC>> bcs(int i) const;

    /// Return trial space
    std::vector<std::shared_ptr<const FunctionSpace>> trial_space() const;
    std::shared_ptr<const FunctionSpace> trial_space(int i) const;

    /// Return test space
    std::vector<std::shared_ptr<const FunctionSpace>> test_space() const;
    std::shared_ptr<const FunctionSpace> test_space(int i) const;

    /// Return lower bound
    std::vector<std::shared_ptr<const GenericVector>> lower_bound() const;

    /// Return upper bound
    std::vector<std::shared_ptr<const GenericVector>> upper_bound() const;

    /// Check whether Jacobian has been defined
    bool has_jacobian() const;

    /// Check whether lower bound has been defined
    bool has_lower_bound() const;

    /// Check whether upper bound have has defined
    bool has_upper_bound() const;

  private:

    // Check forms
    void check_forms() const;

    // Build the necessary mappings between submeshes
    void build_mappings();

    // The residual form
    form_list_type _residual;

    // The Jacobian form (pointer may be null if not provided)
    form_list_type _jacobian;

    // The solution
    std::vector<std::shared_ptr<Function>> _u;

    // The boundary conditions
    std::vector<std::vector<std::shared_ptr<const DirichletBC>>> _bcs;

    // The lower and upper bounds (pointers may be null if not
    // provided)
    std::vector<std::shared_ptr<const GenericVector>> _lb;
    std::vector<std::shared_ptr<const GenericVector>> _ub;
  };

}
    
#endif
