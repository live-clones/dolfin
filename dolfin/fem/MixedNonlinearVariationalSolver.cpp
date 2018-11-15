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
  std::vector<std::shared_ptr<const IndexMap>> index_maps(2);
  MPI_Comm comm = u[0]->vector()->mpi_comm(); // TO CHECK (not tested in //)

  for (size_t i=0; i<u.size(); ++i)
  {
    dolfin_assert(u[i]->vector());
    dolfin_assert(u[i]->function_space()->mesh());
    us.push_back(u[i]->vector());
    
    index_maps[0] = u[i]->function_space()->dofmap().get()->index_map();
    // Create the layout of the jacobian so that we can have J.init_vectors()
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
  // FIXME : We should not need to create a new PETScNestMatrix here
  PETScNestMatrix J(nonlinear_problem->_Js);
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
    // FIXME - TEMPORARY : convergence_criterion set to incremental
    // Default criterion is residual but need to call converged()
    // and this function is not "compatible" with MixedNLSolver (yet)
    newton_solver->parameters.update(parameters("newton_solver"));
    newton_solver->parameters.remove("convergence_criterion");
    newton_solver->parameters.add("convergence_criterion", "incremental");

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
  std::cout << "[MixedNonlinearDiscreteProblem::F] - Implementation in progress... " << std::endl;
  const auto bform = _problem->residual_form();
  auto u = _problem->solution();
  auto bcs = _problem->bcs();

  std::vector<std::shared_ptr<GenericVector>> bs;
  std::vector<std::shared_ptr<GenericVector>> us;
  
  bool has_ufc_form = false;
  std::vector<std::shared_ptr<const IndexMap>> index_maps(2);
  MPI_Comm comm = u[0]->vector()->mpi_comm(); // TO CHECK (not tested in //)

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
    
    if(!has_ufc_form)
    {
      auto tensor_layout = _b->factory().create_layout(comm, 1);
      tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);
      _b->init(*tensor_layout);
      _b->zero();
    }
    
#if 0 // Create a dedicated function here ?
    // Get boundary values for Dirichlet BC
    std::vector<DirichletBC::Map> boundary_values(bcs[i].size());
    for (unsigned int c = 0; c != bcs[i].size(); ++c)
    {
      bcs[i][c]->get_boundary_values(boundary_values[c]);
      if (MPI::size(comm) > 1 && bcs[i][c]->method() != "pointwise")
	bcs[i][c]->gather(boundary_values[c]);
    }

    // Impose Dirichlet BC for the current block (if needed)
    if(!bcs[i].empty())
    {
      const std::size_t num_bc_dofs = boundary_values[0].size();
      std::vector<dolfin::la_index> bc_indices;
      std::vector<double> bc_values;
      bc_indices.reserve(num_bc_dofs);
      bc_values.reserve(num_bc_dofs);

      // Build list of boundary dofs and values
      for (const auto &bv : boundary_values[0])
      {
	bc_indices.push_back(bv.first);
	bc_values.push_back(bv.second);
      }
    
      // Modify bc values
      std::vector<double> x_values(num_bc_dofs);
      u[i]->vector()->get_local(x_values.data(), num_bc_dofs, bc_indices.data());
      for (std::size_t i = 0; i < num_bc_dofs; i++)
	boundary_values[0][bc_indices[i]] = x_values[i] - bc_values[i];
    }
    // Need to change the matrix and the vecctor here from the boundary values that we got
    // Apply_bc is defined in the SystemAssembler class
    // apply_bc(ufc[0]->macro_A.data(), ufc[1]->macro_A.data(), boundary_values,
    // 	     mdofs0, mdofs1);
    
#endif
    us.push_back(u[i]->vector());
    bs.push_back(_b);
  }
  // Combine the vectors
  // FIXME : We should not need to create a new PETScNestMatrix here
  PETScNestMatrix J(_Js);
  J.init_vectors(b, bs);
  #if 0
  for(size_t i=0; i<b.size(); i++)
    std::cout << "[MixedNLSolver] b[" << i << "] = " << b[i] << std::endl;
  #endif
}
//-----------------------------------------------------------------------------
void MixedNonlinearVariationalSolver::
MixedNonlinearDiscreteProblem::J(GenericMatrix& A, const GenericVector& x)
{
  std::cout << "[MixedNonlinearDiscreteProblem::J] - Implementation in progress... " << std::endl;
  const auto jform = _problem->jacobian_form();
  auto u = _problem->solution();
  auto bcs = _problem->bcs();

  //std::vector<std::shared_ptr<GenericMatrix>> Js;
  std::vector<Mat> petsc_mats(u.size()*u.size());
  
  bool has_ufc_form = false;
  std::vector<std::shared_ptr<const IndexMap>> index_maps(2);
  MPI_Comm comm = u[0]->vector()->mpi_comm();

  for (size_t i=0; i<u.size(); ++i)
  {
    index_maps[0] = u[i]->function_space()->dofmap().get()->index_map();
    for (size_t j=0; j<u.size(); ++j)
    {
      has_ufc_form = false;
      std::shared_ptr<GenericMatrix> _J = u[j]->vector()->factory().create_matrix(comm);
      index_maps[1] = u[j]->function_space()->dofmap().get()->index_map();

      for(size_t k=0; k<jform[i*u.size() + j].size(); ++k)
      {
	if(jform[i*u.size() + j][k]->ufc_form()) // If block(i,j) not empty
	{
	  has_ufc_form = true;
	  assemble_mixed(*(_J), *(jform[i*u.size() + j][k]), bool(k>0));
	}
      }

      if(!has_ufc_form)
      {
	auto tensor_layout = _J->factory().create_layout(comm, 2);
	tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);
	_J->init(*tensor_layout);
	_J->zero();
	_J->apply("add");
      }

      _Js[i*u.size() + j] = _J;
      petsc_mats[i*u.size() + j] = as_type<const PETScMatrix>(_J)->mat();
    }
  }
  // TODO : Impose the eventual Dirichlet BC here !

  // Initialize or update the nested matrix
  as_type<PETScMatrix>(A).set_nest(petsc_mats);
  
# if 0 // /!\ DEBUG /!\ Print each blocks to a file
  Mat A00,A01,A10,A11;
  std::vector<IS> is(2);
  MatNestGetISs(as_type<const PETScMatrix>(A).mat(), is.data(), NULL);
  MatGetLocalSubMatrix(as_type<const PETScMatrix>(A).mat(),is[0],is[0],&A00);
  MatGetLocalSubMatrix(as_type<const PETScMatrix>(A).mat(),is[0],is[1],&A01);
  MatGetLocalSubMatrix(as_type<const PETScMatrix>(A).mat(),is[1],is[0],&A10);
  MatGetLocalSubMatrix(as_type<const PETScMatrix>(A).mat(),is[1],is[1],&A11);

  PetscViewer viewer;
  //PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerASCIIOpen(MPI_COMM_WORLD, "A00.m", &viewer);
  MatView(A00,viewer);
  PetscViewerASCIIOpen(MPI_COMM_WORLD, "A01.m", &viewer);
  MatView(A01,viewer);
  PetscViewerASCIIOpen(MPI_COMM_WORLD, "A10.m", &viewer);
  MatView(A10,viewer);
  PetscViewerASCIIOpen(MPI_COMM_WORLD, "A11.m", &viewer);
  MatView(A11,viewer);
#endif
}
