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
// Modified by Cecile Daversin-Catty, 2017.
//
// First added:  2017-07-21
// Last changed: 2017-07-21

#ifndef __MIXED_LINEAR_VARIATIONAL_PROBLEM_H
#define __MIXED_LINEAR_VARIATIONAL_PROBLEM_H

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

  /// This class represents a mixed linear variational problem:
  ///
  /// Find u = (u_1, ..., u_n) in V = V1 x ... x Vn such that
  ///
  ///     a(u, v) = L(v)  for all v in V^,
  ///
  /// where V is the trial space and V^ is the test space.

  class MixedLinearVariationalProblem : public Hierarchical<MixedLinearVariationalProblem>
  {
  public:
    /// Create mixed linear variational problem
    /// with a list of boundary conditions
#if 0
    MixedLinearVariationalProblem(std::vector<std::shared_ptr<const Form>> a,
				  std::vector<std::shared_ptr<const Form>> L,
				  std::vector<std::shared_ptr<Function>> u,
				  std::vector<std::shared_ptr<const DirichletBC>> bcs);
#else
    typedef std::vector<std::vector<std::shared_ptr<const Form>>> form_list_type;
    MixedLinearVariationalProblem(form_list_type a,
				  form_list_type L,
				  std::vector<std::shared_ptr<Function>> u,
				  std::vector<std::shared_ptr<const DirichletBC>> bcs);
#endif

    /// Return bilinear form
    form_list_type bilinear_form() const;
    std::shared_ptr<const Form> bilinear_form(int i, int j=0) const;

    /// Return linear form
    form_list_type linear_form() const;
    std::shared_ptr<const Form> linear_form(int i, int j=0) const;

    /// Return solution variable
    std::vector<std::shared_ptr<Function>> solution();
    std::shared_ptr<Function> solution(int i);

    /// Return boundary conditions
    std::vector<std::vector<std::shared_ptr<const DirichletBC>>> bcs() const;
    std::vector<std::shared_ptr<const DirichletBC>> bcs(int i) const;

    /// Return trial space
    std::vector<std::shared_ptr<const FunctionSpace>> trial_space() const;
    std::shared_ptr<const FunctionSpace> trial_space(int i) const;

    /// Return test space
    std::vector<std::shared_ptr<const FunctionSpace>> test_space() const;
    std::shared_ptr<const FunctionSpace> test_space(int i) const;

  private:

    // Check forms
    void check_forms() const;

    // The bilinear forms
    form_list_type _a;

    // The linear forms
    form_list_type _l;

    // The solution
    std::vector<std::shared_ptr<Function>> _u;

    // The Dirichlet boundary conditions
    std::vector<std::vector<std::shared_ptr<const DirichletBC>>> _bcs;

  };

}

#endif
