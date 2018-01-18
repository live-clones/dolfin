
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

#include <dolfin/common/NoDeleter.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericLinearAlgebraFactory.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/PETScNestMatrix.h>
#include <dolfin/la/PETScKrylovSolver.h>
#include "MixedAssembler.h"
#include "SystemAssembler.h"
#include "assemble.h"
#include "DirichletBC.h"
#include "Form.h"
#include "MixedLinearVariationalProblem.h"
#include "MixedLinearVariationalSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
MixedLinearVariationalSolver::
MixedLinearVariationalSolver(std::shared_ptr<MixedLinearVariationalProblem> problem)
  : _problem(problem)
{
  // Set parameters
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
void MixedLinearVariationalSolver::solve()
{
  begin("Solving mixed linear variational problem.");

  // Get parameters
  std::string solver_type   = parameters["linear_solver"];
  const std::string pc_type = parameters["preconditioner"];
  // const bool print_rhs      = parameters["print_rhs"];
  const bool symmetric      = parameters["symmetric"];
  // const bool print_matrix   = parameters["print_matrix"];

  // Get problem data
  dolfin_assert(_problem);
  const auto a = _problem->bilinear_form();
  const auto L = _problem->linear_form();
  auto u = _problem->solution();
  auto bcs = _problem->bcs();

  // Create matrices and vectors (for each block)
  std::vector<std::shared_ptr<GenericMatrix>> As;
  std::vector<std::shared_ptr<GenericVector>> bs;
  std::vector<std::shared_ptr<GenericVector>> us;
  for (size_t i=0; i<u.size(); ++i)
  {
    for(size_t j=0; j<L[i].size(); ++j)
      dolfin_assert(L[i][j]);
    dolfin_assert(u[i]);

    dolfin_assert(u[i]->vector());
    us.push_back(u[i]->vector());

    // Create rhs vectors
    MPI_Comm comm = u[0]->vector()->mpi_comm(); // TO CHECK (not tested in //)
    std::shared_ptr<GenericVector> b = u[i]->vector()->factory().create_vector(comm);
    bs.push_back(b);

    // Create lhs matrices
    for (size_t j=0; j<u.size(); ++j)
    {
      bool has_ufc_form = false;
      for(int k=0; k<a[i*u.size() + j].size(); ++k)
      {
	dolfin_assert(a[i*u.size() + j][k]);
	if(a[i*u.size() + j][k]->ufc_form())
	  has_ufc_form = true;
      }

      if(has_ufc_form)
      {
	std::shared_ptr<GenericMatrix> A = u[j]->vector()->factory().create_matrix(comm);
	As.push_back(A);
      }
      else
      {
	As.push_back(NULL);
      }
    }
  }

  // Different assembly depending on whether or not the system is symmetric
  if (symmetric)
  {
    // Check that rhs (L) is not empty
    for (size_t i=0; i<L.size(); ++i)
    {
      for (int j=0; j<L[i].size(); ++j)
      {
	if (!L[i][j]->ufc_form())
	{
	  dolfin_error("MixedLinearVariationalSolver.cpp",
		       "symmetric assembly in linear variational solver",
		       "Empty linear forms cannot be used with symmetric assembly");
	}
      }
    }

    // Assemble linear system and apply boundary conditions
    // FIXME : Update the SystemAssembler
    std::cout << "NOT YET IMPLEMENTED" << std::endl;
#if 0
    SystemAssembler assembler(a, L, bcs);
    assembler.assemble(As,bs);
#endif
  }
  else
  {
    // Assemble linear system
    for (size_t i=0; i<u.size(); ++i)
    {
      for (size_t j=0; j<u.size(); ++j)
      {
	// Block-by-block assembly
	for(int k=0; k<a[i*u.size() + j].size(); ++k)
	{
	  if(a[i*u.size() + j][k]->ufc_form()) // If block(i,j) not empty
	    assemble_mixed(*(As[i*u.size() + j]), *(a[i*u.size() + j][k]));
	}
      }
      for (int j=0; j<L[i].size(); ++j)
      {
	if (L[i][j]->ufc_form())
	{
	  // Block-by-block assembly
	  assemble_mixed(*(bs[i]), *(L[i][j]));
	}
	else
	{
	  if (L[i][j]->num_coefficients() != 0)
	  {
	    dolfin_error("MixedLinearVariationalSolver.cpp",
			 "assemble linear form in linear variational solver",
			 "Empty linear forms cannot have coefficient");
	  }
	  for (int k=0; k<u.size(); ++k)
	  {
	    bool has_ufc_form = false;
	    for(int l=0; l<a[i*u.size() + k].size(); ++l)
	    {
	      if(a[i*u.size() + k][l]->ufc_form())
		has_ufc_form = true;
	    }
	    if(has_ufc_form) // If block(i,j) not empty
	      As[i*u.size() + k]->init_vector(*(bs[i]), 0);
	  }
	}
      }
    }

    // Apply boundary conditions
    for (size_t s=0; s<u.size(); ++s) //Number of blocks/subdomains
    {
      for (std::size_t i = 0; i < bcs[s].size(); i++)
      {
	  dolfin_assert(bcs[s][i]);
	  // Apply condition to diagonal blocks
	  bcs[s][i]->apply(*(As[s*u.size() + s]), *(bs[s]));
      }
    }
  }

  // Combine the matrices
  PETScNestMatrix A(As);

  auto comm = u[0]->vector()->mpi_comm();
  auto x = u[0]->vector()->factory().create_vector(comm);
  auto b = u[0]->vector()->factory().create_vector(comm);
  
  // Combine the vectors (solution and rhs)
  A.init_vectors(*x, us);
  A.init_vectors(*b, bs);

  // Solve linear system
  // NOTE : Direct solvers cannot be used with PETScNestMatrix.
  // (For now,) we force the use of a krylov solver
  // (options given by the user for linear_solver are not taken into account)

  // Get list of available preconditioners
  std::map<std::string, std::string>
    preconditioners = u[0]->vector()->factory().krylov_solver_preconditioners();

  // Adjust iterative solver type
  if (symmetric)
    solver_type = "cg";
  else
    solver_type = "gmres";

  if (pc_type != "default"
      && !LinearSolver::in_list(pc_type, preconditioners))
  {
    dolfin_error("MixedLinearVariationalSolver.cpp",
		 "solve mixed linear system",
		 "Unknown preconditioner method \"%s\". "
		 "Use list_krylov_solver_preconditioners() to list available methods",
		 pc_type.c_str());
  }

  PETScKrylovSolver solver(comm, solver_type, pc_type);
  solver.parameters.update(parameters("krylov_solver"));

  solver.set_operator(A);
  solver.solve(*x,*b);

  end();
}
//-----------------------------------------------------------------------------
