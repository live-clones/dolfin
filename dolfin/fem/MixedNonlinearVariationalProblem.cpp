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

#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericVector.h>
#include "Form.h"
#include "DirichletBC.h"
#include "MixedNonlinearVariationalProblem.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
MixedNonlinearVariationalProblem::MixedNonlinearVariationalProblem(
  MixedNonlinearVariationalProblem::form_list_type F,
  std::vector<std::shared_ptr<Function>> u,
  std::vector<std::shared_ptr<const DirichletBC>> bcs,
  MixedNonlinearVariationalProblem::form_list_type J)
  : Hierarchical<MixedNonlinearVariationalProblem>(*this), _residual(F),
  _jacobian(J), _u(u)
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
  if (_residual.size() != _u.size())
      dolfin_error("MixedNonlinearVariationalProblem.cpp",
		   "define mixed nonlinear variational problem F(u, v) == 0 for all v",
		   "Number of blocks in F and solution are inconsistent");
  
  // Check forms
  check_forms();

  // Build the necessary mappings
  build_mappings();
}
//-----------------------------------------------------------------------------
MixedNonlinearVariationalProblem::form_list_type
MixedNonlinearVariationalProblem::residual_form() const
{
  return _residual;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Form>
MixedNonlinearVariationalProblem::residual_form(int i, int j) const
{
  return _residual[i][j];
}
//-----------------------------------------------------------------------------
MixedNonlinearVariationalProblem::form_list_type
MixedNonlinearVariationalProblem::jacobian_form() const
{
  return _jacobian;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Form>
MixedNonlinearVariationalProblem::jacobian_form(int i, int j) const
{
  return _jacobian[i][j];
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<Function>>
MixedNonlinearVariationalProblem::solution()
{
  return _u;
}
//-----------------------------------------------------------------------------
const std::vector<std::shared_ptr<Function>>
MixedNonlinearVariationalProblem::solution() const
{
  return _u;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Function>
MixedNonlinearVariationalProblem::solution(int i)
{
  return _u[i];
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Function>
MixedNonlinearVariationalProblem::solution(int i) const
{
  return _u[i];
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::shared_ptr<const DirichletBC>>>
MixedNonlinearVariationalProblem::bcs() const
{
  return _bcs;
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const DirichletBC>>
MixedNonlinearVariationalProblem::bcs(int i) const
{
  return _bcs[i];
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const FunctionSpace>>
MixedNonlinearVariationalProblem::trial_space() const
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
MixedNonlinearVariationalProblem::test_space() const
{
  std::vector<std::shared_ptr<const FunctionSpace>> test_spaces;
  for (size_t i=0; i<_residual.size(); ++i)
  {
    for (size_t j=0; j<_residual[i].size(); ++j)
      dolfin_assert(_residual[i][j]);
    test_spaces.push_back(_residual[i][0]->function_space(0));
  }
  return test_spaces;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace>
MixedNonlinearVariationalProblem::trial_space(int i) const
{
  dolfin_assert(_u[i]);
  return _u[i]->function_space();
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace>
MixedNonlinearVariationalProblem::test_space(int i) const
{
  for (size_t j=0; j<_residual[i].size(); ++j)
    dolfin_assert(_residual[i][j]);
  return _residual[i][0]->function_space(0);
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const GenericVector>>
MixedNonlinearVariationalProblem::lower_bound() const
{
  for(size_t i=0; i<_lb.size(); ++i)
    dolfin_assert(_lb[i]);
  return _lb;
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const GenericVector>>
MixedNonlinearVariationalProblem::upper_bound() const
{
  for(size_t i=0; i<_ub.size(); ++i)
    dolfin_assert(_ub[i]);
  return _ub;
}
//-----------------------------------------------------------------------------
bool MixedNonlinearVariationalProblem::has_jacobian() const
{
  return !_jacobian.empty();
}
//-----------------------------------------------------------------------------
bool MixedNonlinearVariationalProblem::has_lower_bound() const
{
  return !_lb.empty();
}
//-----------------------------------------------------------------------------
bool MixedNonlinearVariationalProblem::has_upper_bound() const
{
  return !_ub.empty();
}
//-----------------------------------------------------------------------------
void
MixedNonlinearVariationalProblem::build_mappings()
{
  // Residual
  for (size_t i=0; i<_residual.size(); ++i)
  {
    for(size_t j=0; j<_residual[i].size(); ++j)
    {
      // Get integration mesh topology
      if(!_residual[i][j]->ufc_form())
	break;

      // Get meshes associated with the test function
      auto mesh0 = _residual[i][j]->function_space(0)->mesh();
      if(_residual[i][j]->mesh()->id() != mesh0->id()
         && _residual[i][j]->mesh()->topology().mapping().count(mesh0->id()) == 0)
      {
	std::cout << "Build mapping between " << _residual[i][j]->mesh()->id()
		  << " and " << mesh0->id() << std::endl;
        try
        {
          _residual[i][j]->mesh()->build_mapping(mesh0);
        }
        catch (const std::exception& e)
        {
          // If no common parent then restrict the integration domain
          std::cout << e.what() << ", restricting integration domain" << std::endl;
          _residual[i][j]->mesh() = mesh0;
        }
      }

      // Get meshes associated with coefficients (trial)
      auto coefficients = _residual[i][j]->coefficients();
      for(size_t k=0; k<coefficients.size(); ++k)
      {
        if(coefficients[k]->function_space())
        {
          auto mesh1  = coefficients[k]->function_space()->mesh();
          if(_residual[i][j]->mesh()->id() != mesh1->id()
             && _residual[i][j]->mesh()->topology().mapping().count(mesh1->id()) == 0)
          {
            std::cout << "Build mapping between " << _residual[i][j]->mesh()->id()
                      << " and " << mesh1->id() << std::endl;
            try
            {
              _residual[i][j]->mesh()->build_mapping(mesh1);
            }
            catch (const std::exception& e)
            {
              // If no common parent then restrict the integration domain
              std::cout << e.what() << ", restricting integration domain" << std::endl;
              _residual[i][j]->mesh() = mesh1;
            }
          }
        }
      }
    }
  }

  // Jacobian
  for (size_t i=0; i<_jacobian.size(); ++i)
  {
    for(size_t j=0; j<_jacobian[i].size(); ++j)
    {
      // Get integration mesh topology
      if(!_jacobian[i][j]->ufc_form())
	break;

      auto mesh_mapping = _jacobian[i][j]->mesh()->topology().mapping();
      // Get meshes associated with the test function
      auto mesh0 = _jacobian[i][j]->function_space(0)->mesh();
      if(_jacobian[i][j]->mesh()->id() != mesh0->id() && mesh_mapping.count(mesh0->id()) == 0)
      {
	std::cout << "Build mapping between " << _jacobian[i][j]->mesh()->id()
		  << " and " << mesh0->id() << std::endl;
        try
        {
          _jacobian[i][j]->mesh()->build_mapping(mesh0);
        }
        catch (const std::exception& e)
        {
          // If no common parent then restrict the integration domain
          std::cout << e.what() << ", restricting integration domain" << std::endl;
          _jacobian[i][j]->mesh() = mesh0;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------

void
MixedNonlinearVariationalProblem::check_forms() const
{
  std::cout << "[MixedNonlinearVariationalProblem]::check_forms - To be implemented" << std::endl;
}
