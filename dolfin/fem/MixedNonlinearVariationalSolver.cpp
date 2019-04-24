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
// First added:  2018-10-03
// Last changed: 2018-10-03

#include <dolfin/common/NoDeleter.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericLinearAlgebraFactory.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/PETScNestMatrix.h>
#include <dolfin/la/PETScKrylovSolver.h>
#include <dolfin/la/PETScPreconditioner.h>
#include "MixedAssembler.h"
#include "SystemAssembler.h"
#include "assemble.h"
#include "DirichletBC.h"
#include "Form.h"
#include "MixedNonlinearVariationalProblem.h"
#include "MixedNonlinearVariationalSolver.h"
#include "GenericDofMap.h"

#include <dolfin/la/PETScOptions.h>
using namespace dolfin;

MixedNonlinearVariationalSolver::
MixedNonlinearVariationalSolver()
{
  // Set parameters
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
MixedNonlinearVariationalSolver::
MixedNonlinearVariationalSolver(std::shared_ptr<MixedNonlinearVariationalProblem> problem)
  : _problem(problem)
{
  // Set parameters
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool> MixedNonlinearVariationalSolver::solve()
{
  begin("Solving mixed nonlinear variational problem.");

  // Check that the Jacobian has been defined
  dolfin_assert(_problem);
  if (!_problem->has_jacobian())
  {
    dolfin_error("NonlinearVariationalSolver.cpp",
                 "solve nonlinear variational problem",
                 "The Jacobian form has not been defined");
  }
  // Check that the dolfin is configured with petsc is bounds are set
#ifndef HAS_PETSC
  if (_problem->has_lower_bound() || _problem->has_upper_bound())
  {
    dolfin_error("NonlinearVariationalSolver.cpp",
                 "solve nonlinear variational problem",
                 "Needs PETSc to solve bound constrained problems");
  }
#endif

  // Get problem data
  dolfin_assert(_problem);
  auto u = _problem->solution();

  // Create discrete nonlinear problem - assembly
  if (!nonlinear_problem)
  {
    nonlinear_problem = std::make_shared<MixedNonlinearDiscreteProblem>(_problem,
    									reference_to_no_delete_pointer(*this));
  } 

  std::vector<std::shared_ptr<GenericVector>> us;
  MPI_Comm comm = u[0]->vector()->mpi_comm();
  std::vector<std::shared_ptr<const IndexMap>> index_maps(2);

  for (size_t i=0; i<u.size(); ++i)
  {
    dolfin_assert(u[i]->vector());
    us.push_back(u[i]->vector());

    dolfin_assert(u[i]->function_space()->mesh());
    index_maps[0] = u[i]->function_space()->dofmap().get()->index_map();

    // Initialize _bs[i] with the appropriate layout for each block i
    std::shared_ptr<GenericVector> _b = u[i]->vector()->factory().create_vector(comm);
    auto tensor_layout = _b->factory().create_layout(comm, 1);
    tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);
    _b->init(*tensor_layout);
    _b->zero();
    nonlinear_problem->_bs.push_back(_b);

    // Initialize _Js[i,j] with the appropriate layout for each block (i,j)
    for (size_t j=0; j<u.size(); ++j)
    {
      std::shared_ptr<GenericMatrix> _J = u[j]->vector()->factory().create_matrix(comm);
      index_maps[1] = u[j]->function_space()->dofmap().get()->index_map();
      auto tensor_layout = _J->factory().create_layout(comm, 2);
      tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);
      _J->init(*tensor_layout);
      _J->zero();
      _J->apply("add");
      nonlinear_problem->_Js.push_back(_J);
    }
  }

  auto x = u[0]->vector()->factory().create_vector(comm);
  PETScNestMatrix J(nonlinear_problem->_Js);
  // Recombine x = [ u[i]->vector for each block i ]
  // for(auto &u : us) //TO CHECK
  //   as_type<PETScVector>(u)->update_ghost_values();
  J.init_vectors(*x, us);

  std::pair<std::size_t, bool> ret;
  if (std::string(parameters["nonlinear_solver"]) == "newton")
  {
    if (_problem->has_lower_bound() && _problem->has_upper_bound())
    {
      dolfin_error("MixedNonlinearVariationalSolver.cpp",
                   "solve nonlinear variational problem",
                   "Set the \"nonlinear_solver\" parameter to \"snes\" or remove bounds");
    }

    // Create Newton solver and set parameters
    MPI_Comm comm = u[0]->function_space()->mesh()->mpi_comm();
    if (!newton_solver)
    {
      newton_solver = std::make_shared<NewtonSolver>(comm);
    }

    // Pass parameters to Newton solver
    dolfin_assert(newton_solver);
    newton_solver->parameters.update(parameters("newton_solver"));

    // Solve nonlinear problem using Newton's method    
    dolfin_assert(nonlinear_problem);
    ret = newton_solver->solve(*nonlinear_problem, *x);
  }
#ifdef HAS_PETSC
  // FIXME : SNES Solver need to be setup (SNES norm is zero at 1st iteration) 
  else if (std::string(parameters["nonlinear_solver"]) == "snes")
  {
    dolfin_error("MixedNonlinearVariationalSolver.cpp",
                 "solve nonlinear variational problem",
                 "SNES solver is not (yet) available for mixed problems (use newton)");
    
    // Create SNES solver and set parameters
    MPI_Comm comm = u[0]->function_space()->mesh()->mpi_comm();
    if (!snes_solver)
    {
      snes_solver = std::make_shared<PETScSNESSolver>(comm);
    }
    snes_solver->parameters.update(parameters("snes_solver"));

    // Solve nonlinear problem using PETSc's SNES
    dolfin_assert(nonlinear_problem);
    if (_problem->has_lower_bound() && _problem->has_upper_bound())
    {
      ret = snes_solver->solve(*nonlinear_problem, *x,
                               *_problem->lower_bound()[0],
                               *_problem->upper_bound()[0]);
    }
    else
    {
      ret = snes_solver->solve(*nonlinear_problem, *x);
    }
  }
#endif
  else
  {
    dolfin_error("MixedNonlinearVariationalSolver.cpp",
                 "solve nonlinear variational problem",
                 "Unknown nonlinear solver type");
  }

  end();
  return ret;
}
//-----------------------------------------------------------------------------
// Implementation of NonlinearDiscreteProblem
//-----------------------------------------------------------------------------
MixedNonlinearVariationalSolver::MixedNonlinearDiscreteProblem::
MixedNonlinearDiscreteProblem(std::shared_ptr<const MixedNonlinearVariationalProblem> problem,
			      std::shared_ptr<const MixedNonlinearVariationalSolver> solver)
  : _problem(problem), _solver(solver)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MixedNonlinearVariationalSolver::MixedNonlinearDiscreteProblem::~MixedNonlinearDiscreteProblem()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MixedNonlinearVariationalSolver::
MixedNonlinearDiscreteProblem::F(GenericVector& b, const GenericVector& x)
{
  const auto bform = _problem->residual_form();
  auto u = _problem->solution();
  auto bcs = _problem->bcs();
  
  bool has_ufc_form = false;
  MPI_Comm comm = u[0]->vector()->mpi_comm();

  for(size_t i=0; i<u.size(); ++i)
  {
    has_ufc_form = false;
    std::shared_ptr<GenericVector> _b = u[i]->vector()->factory().create_vector(comm);

    for(size_t j=0; j<bform[i].size(); ++j)
    {
      if (bform[i][j]->ufc_form())
      {
	has_ufc_form = true;
	assemble_mixed(*(_b), *(bform[i][j]), bool(j>0));
      }
    }

    // Update _bs
    if(has_ufc_form)
      _bs[i] = _b;
  }

  // Boundary conditions
  PETScNestMatrix J(_Js);
  std::vector<IS> is(u.size());
  MatNestGetISs(J.mat(), is.data(), NULL);

  for(size_t i=0; i<u.size(); ++i)
  {
    // Extract the subvector x[i] of x corresponding to block i
    Vec __x;
    VecGetSubVector(as_type<const PETScVector>(x).vec(),is[i], &__x);
    PETScVector _x(__x);
    
    for (size_t c = 0; c < bcs[i].size(); c++)
    {
      dolfin_assert(bcs[i][c]);
      bcs[i][c]->apply(*_bs[i], _x);
    }
  }

  // Update b from _bs[]
  J.init_vectors(b, _bs);
}
//-----------------------------------------------------------------------------
void MixedNonlinearVariationalSolver::
MixedNonlinearDiscreteProblem::J(GenericMatrix& A, const GenericVector& x)
{
  const auto jform = _problem->jacobian_form();
  auto u = _problem->solution();
  auto bcs = _problem->bcs();

  MPI_Comm comm = u[0]->vector()->mpi_comm();

  for (size_t i=0; i<u.size(); ++i)
  {
    for (size_t j=0; j<u.size(); ++j)
    {
      for(size_t k=0; k<jform[i*u.size() + j].size(); ++k)
      {
	if(jform[i*u.size() + j][k]->ufc_form()) // If block(i,j) not empty
	{
          _Js[i*u.size() + j] = u[j]->vector()->factory().create_matrix(comm);
          assemble_mixed(*(_Js[i*u.size() + j]), *(jform[i*u.size() + j][k]), bool(k>0));
	}
      }
    }
  }

  // Boundary conditions
  PETScNestMatrix J(_Js);
  std::vector<IS> is(u.size());
  MatNestGetISs(J.mat(), is.data(), NULL);

  for(size_t i=0; i<u.size(); ++i)
  {
    Vec __x;
    VecGetSubVector(as_type<const PETScVector>(x).vec(),is[i], &__x);
    PETScVector _x(__x);

    for (size_t c = 0; c < bcs[i].size(); c++)
    {
      dolfin_assert(bcs[i][c]);
      for (size_t j=0; j<u.size(); ++j)
      {
        if(i!=j && jform[i*u.size() + j][0]->ufc_form())
          bcs[i][c]->zero_columns(*(_Js[i*u.size() + j]), *(_bs[i]), 0, true);
        else
          bcs[i][c]->apply(*(_Js[i*u.size() + i]), *(_bs[i]), _x);
      }
    }
  }

  // Update A from _Js, as a PETScNestMatrix
  // FIXME : We should not need petsc_mats and find a way to call set_nest on _Js
  std::vector<Mat> petsc_mats(u.size()*u.size());
  for(size_t i=0; i<u.size(); ++i)
    for(size_t j=0; j<u.size(); ++j)
      petsc_mats[i*u.size() + j] = as_type<PETScMatrix>(_Js[i*u.size() + j])->mat();
  as_type<PETScMatrix>(A).set_nest(petsc_mats);
}
