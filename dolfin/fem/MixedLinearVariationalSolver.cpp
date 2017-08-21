
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
  const bool print_rhs      = parameters["print_rhs"];
  const bool symmetric      = parameters["symmetric"];
  const bool print_matrix   = parameters["print_matrix"];

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
  for (int i=0; i<u.size(); ++i)
  {
    dolfin_assert(L[i]);
    dolfin_assert(u[i]);

    dolfin_assert(u[i]->vector());
    us.push_back(u[i]->vector());

    // Create rhs vectors
    MPI_Comm comm = u[i]->vector()->mpi_comm(); // TO CHECK
    std::shared_ptr<GenericVector> b = u[i]->vector()->factory().create_vector(comm);
    bs.push_back(b);

    // Create lhs matrices
    for (int j=0; j<u.size(); ++j)
    {
      dolfin_assert(a[i*u.size() + j]);
      std::shared_ptr<GenericMatrix> A = u[j]->vector()->factory().create_matrix(comm);
      if(a[i*u.size() + j]->ufc_form())
	As.push_back(A);
      else
	As.push_back(NULL);
    }
  }

  // Different assembly depending on whether or not the system is symmetric
  if (symmetric)
  {
    // Check that rhs (L) is not empty
    for (int i=0; i<L.size(); ++i)
    {
	if (!L[i]->ufc_form())
	{
	    dolfin_error("MixedLinearVariationalSolver.cpp",
			 "symmetric assembly in linear variational solver",
			 "Empty linear forms cannot be used with symmetric assembly");
	}
    }

    // Assemble linear system and apply boundary conditions
    SystemAssembler assembler(a, L, bcs);
    assembler.assemble(As,bs);
  }
  else
  {
    // Assemble linear system
    for (int i=0; i<u.size(); ++i)
    {
      for (int j=0; j<u.size(); ++j)
      {
	// Block-by-block assembly
	// std::cout << "a block (" << i << "," << j << ")"<< std::endl;
	if(a[i*u.size() + j]->ufc_form()) // If block(i,j) not empty
	  assemble_mixed(*(As[i*u.size() + j]), *(a[i*u.size() + j]));
      }
      if (L[i]->ufc_form())
      {
	// Block-by-block assembly
	assemble_mixed(*(bs[i]), *(L[i]));
      }
      else
      {
	if (L[i]->num_coefficients() != 0)
	{
	  dolfin_error("MixedLinearVariationalSolver.cpp",
		       "assemble linear form in linear variational solver",
		       "Empty linear forms cannot have coefficient");
	}
	for (int j=0; j<u.size(); ++j)
	  if(a[i*u.size() + j]->ufc_form()) // If block(i,j) not empty
	    As[i*u.size() + j]->init_vector(*(bs[i]), 0);
      }
    }

    // Apply boundary conditions
    for (int s=0; s<u.size(); ++s) //Number of blocks/subdomains
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
  A.init_vectors(*x, us);
  A.init_vectors(*b, bs);

  // Solve linear system
  // TEMPORARY : Force PETScKrylovSolver (need to be carefully configured)
  PETScKrylovSolver solver(comm, solver_type, pc_type);
  solver.parameters.update(parameters("krylov_solver"));

  solver.set_operators(A,A);
  solver.solve(*x,*b);

  end();
}
//-----------------------------------------------------------------------------
