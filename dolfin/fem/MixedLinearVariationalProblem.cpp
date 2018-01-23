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

#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include "Form.h"
#include "DirichletBC.h"
#include "MixedLinearVariationalProblem.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
#if 0
MixedLinearVariationalProblem::MixedLinearVariationalProblem(
  std::vector<std::shared_ptr<const Form>> a,
  std::vector<std::shared_ptr<const Form>> L,
  std::vector<std::shared_ptr<Function>> u,
  std::vector<std::shared_ptr<const DirichletBC>> bcs)
  : Hierarchical<MixedLinearVariationalProblem>(*this), _a(a), _l(L), _u(u)
#else
MixedLinearVariationalProblem::MixedLinearVariationalProblem(
  MixedLinearVariationalProblem::form_list_type a,
  MixedLinearVariationalProblem::form_list_type L,
  std::vector<std::shared_ptr<Function>> u,
  std::vector<std::shared_ptr<const DirichletBC>> bcs)
  : Hierarchical<MixedLinearVariationalProblem>(*this), _a(a), _l(L), _u(u)  
#endif
{
  // Initialize each sub vectors
  for (size_t i=0; i<u.size(); ++i)
      _bcs.push_back( std::vector<std::shared_ptr<const DirichletBC>>() );
  
  // Sort conditions from the function space
  for (size_t i=0; i<bcs.size(); ++i)
  {
      for (size_t j=0; j<_u.size(); ++j)
      {
	  if (_u[j]->in(*bcs[i]->function_space()))
	      _bcs[j].push_back(bcs[i]);
      }
      
  }

  // Check size
  if (_a.size() != _u.size()*_u.size() || _l.size() != _u.size())
      dolfin_error("MixedLinearVariationalProblem.cpp",
		   "define mixed linear variational problem a(u, v) == L(v) for all v",
		   "Number of blocks in rhs, lhs, solution are inconsistent");
  
  // Check forms
  // FIXME
  // check_forms();
}
//-----------------------------------------------------------------------------
MixedLinearVariationalProblem::form_list_type
MixedLinearVariationalProblem::bilinear_form() const
{
  return _a;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Form>
MixedLinearVariationalProblem::bilinear_form(int i, int j) const
{
  return _a[i][j];
}
//-----------------------------------------------------------------------------
MixedLinearVariationalProblem::form_list_type
MixedLinearVariationalProblem::linear_form() const
{
  return _l;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Form>
MixedLinearVariationalProblem::linear_form(int i, int j) const
{
  return _l[i][j];
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<Function>>
MixedLinearVariationalProblem::solution()
{
  return _u;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Function>
MixedLinearVariationalProblem::solution(int i)
{
  return _u[i];
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::shared_ptr<const DirichletBC>>>
MixedLinearVariationalProblem::bcs() const
{
  return _bcs;
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const DirichletBC>>
MixedLinearVariationalProblem::bcs(int i) const
{
  return _bcs[i];
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const FunctionSpace>>
MixedLinearVariationalProblem::trial_space() const
{
  std::vector<std::shared_ptr<const FunctionSpace>> trial_spaces;
  for (size_t i=0; i<_u.size(); ++i)
  {
    dolfin_assert(_u[i]);
    trial_spaces.push_back(_u[i]->function_space());
  }
  return trial_spaces;
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const FunctionSpace>>
MixedLinearVariationalProblem::test_space() const
{
  std::vector<std::shared_ptr<const FunctionSpace>> test_spaces;
  for (size_t i=0; i<_l.size(); ++i)
  {
    for (size_t j=0; j<_l[i].size(); ++j)
      dolfin_assert(_l[i][j]);
    test_spaces.push_back(_l[i][0]->function_space(0));
    // TODO : Check if all have the same FS
  }
  return test_spaces;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace>
MixedLinearVariationalProblem::trial_space(int i) const
{
  dolfin_assert(_u[i]);
  return _u[i]->function_space();
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace>
MixedLinearVariationalProblem::test_space(int i) const
{
  for (size_t j=0; j<_l[i].size(); ++j)
    dolfin_assert(_l[i][j]);
  // FIXME
  return _l[i][0]->function_space(0);
}
//-----------------------------------------------------------------------------

void
MixedLinearVariationalProblem::check_forms() const
{
  // Check rank of bilinear form a
  for (size_t i=0; i<_a.size(); ++i)
  {
    for(int j=0; j<_a[i].size(); ++j)
    {
      dolfin_assert(_a[i][j]);
      if (_a[i][j]->rank() != 2)
      {
	dolfin_error("MixedLinearVariationalProblem.cpp",
		     "define mixed linear variational problem a(u, v) == L(v) for all v",
		     "Expecting the left-hand side (block %d, domain %d) to be a bilinear form (not rank %d)",
		     i, j, _a[i][j]->rank());
      }
    }
  }

  // Check rank of i-th linear form L
  for (int i=0; i<_l.size(); ++i)
  {
    for(int j=0; j<_l[i].size(); ++j)
    {
      // Check rank of i-th linear form L
      dolfin_assert(_l[i][j]);
      if (_l[i][j]->rank() != 1)
      {
	  dolfin_error("MixedLinearVariationalProblem.cpp",
		       "define mixed linear variational problem a(u, v) = L(v) for all v",
		       "Expecting the right-hand side (block %d, domain %d) to be a linear form (not rank %d)",
		       i,j,_l[i][j]->rank());
      }
    }

    // Check that function space of solution variable matches trial space (for each problem)
    dolfin_assert(_u[i]);
    for (int j=0; j<_l.size(); ++j)
    {
      for(int k=0; k<_a[i + j*_l.size()].size(); ++k)
      {
	const auto trial_space = _a[i + j*_l.size()][k]->function_space(1);
	dolfin_assert(trial_space);
	// trial_space can be null is the block is empty (Form=None)
	if (trial_space != nullptr && !_u[i]->in(*trial_space))
	{
	  dolfin_error("MixedLinearVariationalProblem.cpp",
		       "define mixed linear variational problem a(u, v) = L(v) for all v",
		       "Expecting the solution variable u to be a member of the trial space");
	}
      }
    }

    // Check that function spaces of bcs are contained in trial space
    for (const auto bc: _bcs[i])
    {
      dolfin_assert(bc);
      const auto bc_space = bc->function_space();
      dolfin_assert(bc_space);
      for (int j=0;  j<_l.size(); ++j)
      {
	for(int k=0; k<_a[i + j*_l.size()].size(); ++k)
	{
	  const auto trial_space = _a[i + j*_l.size()][k]->function_space(1);
	  if (trial_space != nullptr && !trial_space->contains(*bc_space))
	  {
	    dolfin_error("MixedLinearVariationalProblem.cpp",
			 "define mixed linear variational problem a(u, v) = L(v) for all v",
			 "Expecting the boundary conditions to live on (a "
			 "subspace of) the trial space");
	  }
	}
      }
    }
  }
}
//-----------------------------------------------------------------------------
