
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
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericLinearAlgebraFactory.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/PETScNestMatrix.h>
#include <dolfin/la/PETScKrylovSolver.h>
#include <dolfin/la/PETScLUSolver.h>
#include <dolfin/la/PETScPreconditioner.h>
#include "MixedAssembler.h"
#include "SystemAssembler.h"
#include "assemble.h"
#include "DirichletBC.h"
#include "Form.h"
#include "MixedLinearVariationalProblem.h"
#include "MixedLinearVariationalSolver.h"
#include "GenericDofMap.h"

#include <dolfin/la/PETScOptions.h>
using namespace dolfin;

//-----------------------------------------------------------------------------
MixedLinearVariationalSolver::
MixedLinearVariationalSolver()
{
  // Set parameters
  parameters = default_parameters();
}


MixedLinearVariationalSolver::
MixedLinearVariationalSolver(std::shared_ptr<MixedLinearVariationalProblem> problem)
  : _problem(problem)
{
  // Set parameters
  parameters = default_parameters();
}

MixedLinearVariationalSolver::assembled_system_type
MixedLinearVariationalSolver::assemble_system()
{
  const bool symmetric = parameters["symmetric"];

  const auto a = _problem->bilinear_form();
  const auto L = _problem->linear_form();
  auto u = _problem->solution();
  auto bcs = _problem->bcs();
  
  // Create matrices and vectors (for each block)
  std::vector<std::shared_ptr<GenericMatrix>> As;
  std::vector<std::shared_ptr<GenericVector>> bs;
  std::vector<std::shared_ptr<GenericVector>> us;

  std::vector<std::shared_ptr<const IndexMap>> index_maps(2);
  bool has_ufc_form = false;

  MPI_Comm comm = u[0]->vector()->mpi_comm();

  for (size_t i=0; i<u.size(); ++i)
  {
    dolfin_assert(u[i]);
    dolfin_assert(u[i]->vector());
    us.push_back(u[i]->vector());

    // Create rhs vectors
    for(size_t j=0; j<L[i].size(); ++j)
    {
      dolfin_assert(L[i][j]);
      if (L[i][j]->ufc_form())
	has_ufc_form = true;
    }

    std::shared_ptr<GenericVector> b = u[i]->vector()->factory().create_vector(comm);
    index_maps[0] = u[i]->function_space()->dofmap().get()->index_map();

    if(!has_ufc_form)
    {
      auto tensor_layout = b->factory().create_layout(comm, 1);
      tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);
      b->init(*tensor_layout);
      b->zero();
    }
    bs.push_back(b);

    // Create lhs matrices
    for (size_t j=0; j<u.size(); ++j)
    {
      has_ufc_form = false;
      for(size_t k=0; k<a[i*u.size() + j].size(); ++k)
      {
	dolfin_assert(a[i*u.size() + j][k]);
	if(a[i*u.size() + j][k]->ufc_form())
	  has_ufc_form = true;
      }

      std::shared_ptr<GenericMatrix> A = u[j]->vector()->factory().create_matrix(comm);
      index_maps[1] = u[j]->function_space()->dofmap().get()->index_map();

      if(!has_ufc_form)
      {
	auto tensor_layout = A->factory().create_layout(comm, 2);
	tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);
	A->init(*tensor_layout);
	A->zero();
	A->apply("add");
      }
      As.push_back(A);
    }
  }

  // Different assembly depending on whether or not the system is symmetric
  if (symmetric)
  {
    // Check that rhs (L) is not empty
    for (size_t i=0; i<L.size(); ++i)
    {
      for (size_t j=0; j<L[i].size(); ++j)
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
  }
  else
  {
    // Assemble linear system
    for (size_t i=0; i<u.size(); ++i)
    {
      for (size_t j=0; j<u.size(); ++j)
      {
	// Block-by-block assembly
	for(size_t k=0; k<a[i*u.size() + j].size(); ++k)
	{
	  if(a[i*u.size() + j][k]->ufc_form()) // If block(i,j) not empty
	    assemble_mixed(*(As[i*u.size() + j]), *(a[i*u.size() + j][k]), bool(k>0));
	}
      }
      for (size_t j=0; j<L[i].size(); ++j)
      {
	if (L[i][j]->ufc_form())
	{
	  // Block-by-block assembly
	  assemble_mixed(*(bs[i]), *(L[i][j]), bool(j>0));
	}
	else
	{
	  if (L[i][j]->num_coefficients() != 0)
	  {
	    dolfin_error("MixedLinearVariationalSolver.cpp",
			 "assemble linear form in linear variational solver",
			 "Empty linear forms cannot have coefficient");
	  }
	  for (size_t k=0; k<u.size(); ++k)
	  {
	    bool has_ufc_form = false;
	    for(size_t l=0; l<a[i*u.size() + k].size(); ++l)
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
    for (size_t i=0; i<u.size(); ++i)
    {
      // Iterate over bcs for i-th diag block
      for (size_t c = 0; c < bcs[i].size(); c++)
      {
	dolfin_assert(bcs[i][c]);
	for (size_t j=0; j<u.size(); ++j)
	{
	  if(i!=j && a[i*u.size() + j][0]->ufc_form()) // Apply condition to off-diag blocks (zero values)
	    bcs[i][c]->zero_columns(*(As[i*u.size() + j]), *(bs[i]), 0, true);
	  else // Apply condition to diag blocks
	    bcs[i][c]->apply(*(As[i*u.size() + i]), *(bs[i]));
	}
      }
    }
  }
  return std::make_tuple(As,bs,us);
}
//-----------------------------------------------------------------------------
void MixedLinearVariationalSolver::solve()
{
  auto assembled_system = this->assemble_system();
  this->solve(assembled_system);
}
//-----------------------------------------------------------------------------
void MixedLinearVariationalSolver::solve(MixedLinearVariationalSolver::assembled_system_type assembled_system)
{
  begin("Solving mixed linear variational problem.");

  // Get parameters
  std::string solver_type   = parameters["linear_solver"];
  const std::string pc_type = parameters["preconditioner"];

  auto As = std::get<0>(assembled_system);
  auto bs = std::get<1>(assembled_system);
  auto us = std::get<2>(assembled_system);

  // Combine the matrices
  PETScNestMatrix A(As);

  MPI_Comm comm = us[0]->mpi_comm();
  std::shared_ptr<GenericVector> x = us[0]->factory().create_vector(comm);
  std::shared_ptr<GenericVector> b = us[0]->factory().create_vector(comm);
  
  // Combine the vectors (solution and rhs)
  A.init_vectors(*x, us);
  A.init_vectors(*b, bs);

  // Get list of available methods
  std::map<std::string, std::string>
    lu_methods = us[0]->factory().lu_solver_methods();
  std::map<std::string, std::string>
    krylov_methods = us[0]->factory().krylov_solver_methods();
  std::map<std::string, std::string>
    preconditioners = us[0]->factory().krylov_solver_preconditioners();

  if (!LinearSolver::in_list(solver_type, krylov_methods)
      && LinearSolver::in_list(solver_type, lu_methods))
  {
    dolfin_error("LinearVariationalSolver.cpp",
		 "solve linear system",
		 "Unknown solver method \"%s\". "
		 "Use list_linear_solver_methods() to list available methods",
		 solver_type.c_str());
  }

  if (pc_type != "default"
      && !LinearSolver::in_list(pc_type, preconditioners))
  {
    dolfin_error("MixedLinearVariationalSolver.cpp",
		 "solve mixed linear system",
		 "Unknown preconditioner method \"%s\". "
		 "Use list_krylov_solver_preconditioners() to list available methods",
		 pc_type.c_str());
  }

  if (solver_type == "direct" || solver_type == "lu")
  {
    solver_type = "default";
    PETScLUSolver solver(comm, solver_type);
    // Convert from MATNEST to AIJ matrix type
    std::cout << "Converting PETScNestMatrix into AIJ (due to direct solver)" << std::endl;
    A.convert_to_aij();
    solver.set_operator(A);
    solver.solve(*x,*b);
  }
  else
  {
    PETScKrylovSolver solver(comm, solver_type, pc_type);
    solver.parameters.update(parameters("krylov_solver"));
    if(pc_type != "default")
    {
      std::cout << "Converting PETScNestMatrix into AIJ (due to "
                << pc_type  << " preconditioner)" << std::endl;
      A.convert_to_aij();
    }
    solver.set_operator(A);
    solver.solve(*x,*b);
  }
  // Ghost values have to be updated in parallel since they are not updated by solve
  for(auto &u : us)
    as_type<PETScVector>(u)->update_ghost_values();

#if 0 // Configure PCFieldSplit (Not working yet)
  PETScOptions::set("pc_type", "fieldsplit");
  PETScOptions::set("pc_fieldsplit_type", "additive");
  solver.set_from_options();

  std::vector<std::vector<dolfin::la_index>> fields(u.size());
  std::vector<std::string> split_names(u.size());
  for(size_t i=0; i<u.size(); ++i)
  {
    //split_names[i] = "u" + std::to_string(i);
    split_names[i] = std::to_string(i);
    A.get_block_dofs(fields[i], i);

    // Ksp_type
    std::string ksp_type_name = "fieldsplit_" + split_names[i] + "_ksp_type";
    std::string pc_type_name = "fieldsplit_" + split_names[i] + "_pc_type";
    PETScOptions::set(ksp_type_name, solver_type);
    PETScOptions::set(pc_type_name, pc_type);
  }
  //solver.set_from_options();
  PETScPreconditioner::set_fieldsplit(solver, fields, split_names);
  //KSPView(solver.ksp(), PETSC_VIEWER_STDOUT_WORLD);
#endif
  end();
}
