// Copyright (C) 2016 Garth N. Wells and Chris N. Richardson
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
// This demo program solves Stokes's equations
//
//     mu*div(grad(u)) - grad(p) = f
//     div(u) = 0
//
// on a unit cube, using the MatNest facility of PETSc. By doing this, it is possible
// to save memory when creating preconditioners, an essential requirement for computation
// at large scale (over 100 million degrees of freedom).

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include <dolfin.h>

#include<petscsys.h>
#include<petscksp.h>

#include "Stokes.h"

using namespace dolfin;

// Function to compute the near nullspace A00 operator
VectorSpaceBasis build_nullspace_nested(const FunctionSpace& V,
                                                const GenericVector& x)
{
  // Get subspaces
  auto V0 = V.sub(0);
  auto V1 = V.sub(1);
  auto V2 = V.sub(2);

  // Create vectors for nullspace basis
  std::vector<std::shared_ptr<GenericVector>> basis(3);
  for (std::size_t i = 0; i < basis.size(); ++i)
    basis[i] = x.copy();

  // x0, x1, x2 translations
  V0->dofmap()->set(*basis[0], 1.0);
  V1->dofmap()->set(*basis[1], 1.0);
  V2->dofmap()->set(*basis[2], 1.0);

  // Rotations
  /*
  V0->set_x(*basis[3], -1.0, 1);
  V1->set_x(*basis[3],  1.0, 0);

  V0->set_x(*basis[4],  1.0, 2);
  V2->set_x(*basis[4], -1.0, 0);

  V2->set_x(*basis[5],  1.0, 1);
  V1->set_x(*basis[5], -1.0, 2);
  */

  // Apply
  for (std::size_t i = 0; i < basis.size(); ++i)
    basis[i]->apply("add");

  // Create vector space and orthonormalize
  VectorSpaceBasis vector_space(basis);
  vector_space.orthonormalize();
  return vector_space;
}


std::array<std::shared_ptr<Function>, 2>
  solve_nested(std::shared_ptr<Mesh> mesh,
               std::shared_ptr<SubDomain> top_bottom,
               std::shared_ptr<SubDomain> left_edge,
               std::shared_ptr<GenericFunction> flow_velocity,
               std::shared_ptr<GenericFunction> pleft,
               std::shared_ptr<GenericFunction> f,
               std::shared_ptr<GenericFunction> zero,
               bool pc_operator)
{
}


// Function for no-slip boundary condition for velocity
class Source : public Expression
{
public:

  Source() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = x[1];
    values[1] = x[1]*x[2];
    values[2] = 0.0;
  }
};

// Sub domain for top and bottom
class TopBottom : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary;
    // return std::abs(1.0 - x[1]) < DOLFIN_EPS || std::abs(x[1]) < DOLFIN_EPS;
  }
};

class LeftEdge : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS;
  }
};

int main(int argc, char *argv[])
{
  // Application parameters
  Parameters application_parameters("application_parameters");
  application_parameters.add("size", 12);

  // Parse application parameters
  application_parameters.parse(argc, argv);
  const int size = application_parameters["size"];

// Parse PETSc parameters
  parameters.parse(argc, argv);

  PETScOptions::set("ksp_monitor_true_residual");

  // Create mesh
  auto mesh = std::make_shared<UnitCubeMesh>(size, size, size);

  // Velocity BC
  auto top_bottom = std::make_shared<TopBottom>();
  auto flow_velocity = std::make_shared<Constant>(0.0, 0.0, 0.0) ;

  auto f = std::make_shared<Source>();
  auto zero = std::make_shared<Constant>(0.0);

  // Common fieldsplit options
  PETScOptions::set("ksp_type", "minres");
  PETScOptions::set("pc_type", "fieldsplit");
  PETScOptions::set("pc_fieldsplit_type", "additive");

  PETScOptions::set("ksp_view");

  // AMG
  PETScOptions::set("fieldsplit_0_ksp_type", "preonly");

  // "gamg" or "hypre"
  std::string mg = "hypre";

  if (mg == "gamg")
  {
    PETScOptions::set("fieldsplit_0_pc_type", "gamg");
    PETScOptions::set("fieldsplit_0_pc_gamg_coarse_eq_limit", 10000);

    PETScOptions::set("fieldsplit_0_pc_gamg_threshold", 0.05);
    //PETScOptions::set("fieldsplit_0_pc_gamg_square_graph", 2);
    PETScOptions::set("fieldsplit_0_mg_levels_ksp_type", "chebyshev");
    PETScOptions::set("fieldsplit_0_mg_levels_pc_type", "sor");
    PETScOptions::set("fieldsplit_0_mg_levels_ksp_max_it", 3);

    PETScOptions::set("fieldsplit_0_mg_levels_esteig_ksp_type", "cg");
    PETScOptions::set("fieldsplit_0_mg_levels_ksp_chebyshev_esteig_steps", 20);
    PETScOptions::set("fieldsplit_0_mg_levels_ksp_chebyshev_esteig_random");

    PETScOptions::set("fieldsplit_0_mg_coarse_ksp_type", "preonly");
    PETScOptions::set("fieldsplit_0_mg_coarse_pc_type", "lu");
    PETScOptions::set("fieldsplit_0_mg_coarse_pc_factor_mat_solver_package", "mumps");
  }
  else if (mg == "hypre")
  {
    PETScOptions::set("fieldsplit_0_pc_type", "hypre");
    PETScOptions::set("fieldsplit_0_hypre_type", "boomeramg");
    PETScOptions::set("fieldsplit_1_pc_hypre_boomeramg_strong_threshold", 0.5);
  }
  else
    error("Unknown A00 pc");

  PETScOptions::set("fieldsplit_1_ksp_type", "preonly");
  PETScOptions::set("fieldsplit_1_pc_type", "hypre");
  PETScOptions::set("fieldsplit_1_hypre_type", "boomeramg");
  PETScOptions::set("fieldsplit_1_pc_hypre_boomeramg_strong_threshold", 0.5);

  // Create function spaces
  auto V = std::make_shared<Stokes::Form_a00::TestSpace>(mesh);
  auto Q = std::make_shared<Stokes::Form_a11::TestSpace>(mesh);

  // Define variational problem
  auto a00 = std::make_shared<Stokes::Form_a00>(V, V);
  auto L0 = std::make_shared<Stokes::Form_L0>(V);
  L0->f = f;

  auto a01 = std::make_shared<Stokes::Form_a01>(Q, V);
  auto a10 = std::make_shared<Stokes::Form_a10>(V, Q);
  auto a11 = std::make_shared<Stokes::Form_a11>(Q, Q);
  a11->c = zero;

  auto L1 = std::make_shared<Stokes::Form_L1>(Q);
  L1->g = zero;

  // Velocity BC
  auto bcV = std::make_shared<DirichletBC>(V, flow_velocity,
                                                   top_bottom);

  // Assemble all blocks with BCs
  SystemAssembler assemblerA0({a00, a01, a10, a11}, {L0, L1},
                                      {{bcV}, {}});
  auto A00 = std::make_shared<PETScMatrix>();
  auto A01 = std::make_shared<PETScMatrix>();
  auto b0 = std::make_shared<PETScVector>();
  auto A10 = std::make_shared<PETScMatrix>();
  auto b1 = std::make_shared<PETScVector>();

  auto A11 = std::make_shared<PETScMatrix>();
  std::cout << "Assemble A00, A01, A10, A11\n";
  assemblerA0.assemble({A00, A01, A10, A11}, {b0, b1});

  // Assemble a mass matrix (QxQ) for the preconditioner
  auto one = std::make_shared<Constant>(1.0);
  a11->c = one;
  SystemAssembler assemblerP(a11, L1, {});
  auto P11 = std::make_shared<PETScMatrix>();
  assemblerP.assemble(*P11);

  std::cout << "A00:" << A00->size(0) << "x" << A00->size(1)
            << " : " << std::setprecision(10)
            << std::pow(A00->norm("frobenius"), 2) << "\n";
  std::cout << "A01:" << A01->size(0) << "x" << A01->size(1)
            << " : " << std::setprecision(10)
            << std::pow(A01->norm("frobenius"), 2) << "\n";
  std::cout << "A10:" << A10->size(0) << "x" << A10->size(1)
            << " : " << std::setprecision(10)
            << std::pow(A10->norm("frobenius"), 2) << "\n";
  std::cout << "P11:" << P11->size(0) << "x" << P11->size(1)
            << " : " << std::setprecision(10)
            << std::pow(P11->norm("frobenius"), 2) << "\n";

  std::cout << "b0:" << b0->local_range().second - b0->local_range().first
            << " : " << std::setprecision(10)
            << std::pow(b0->norm("l2"), 2) << "\n";
  std::cout << "b1:" << b1->local_range().second - b1->local_range().first
            << " : " << std::setprecision(10)
            << std::pow(b1->norm("l2"), 2) << "\n";

  // Create nested A and P

  // NOTE: We need an empty bottom-right block when using PETSc Schur
  // complement preconditioning.
  std::vector<std::shared_ptr<const GenericMatrix>> Amats, Pmats;
  std::shared_ptr<PETScNestMatrix> Anest, Pnest;

  Amats = {A00, A01, A10, NULL};
  Anest = std::make_shared<PETScNestMatrix>(Amats);

  Pmats = {A00, NULL, NULL, P11};
  Pnest = std::make_shared<PETScNestMatrix>(Pmats);

  // Create VecNest vectors for RHS and solution vector
  auto u0 = std::make_shared<Function>(V);
  auto p0 = std::make_shared<Function>(Q);
  PETScVector xvec, bvec;
  Anest->init_vectors(xvec, {u0->vector(), p0->vector()});
  Anest->init_vectors(bvec, {b0, b1});

  // ---------
  // Create pressure nullspace vector
  auto u0tmp = std::make_shared<Function>(V);
  auto p0tmp = std::make_shared<Function>(Q);
  p0tmp->interpolate(*one);

  Constant one_vec(0.0, 0.0, 0.0);
  u0tmp->interpolate(one_vec);

  PETScVector xvec_nullspace;
  Anest->init_vectors(xvec_nullspace, {u0tmp->vector(), p0tmp->vector()});
  xvec_nullspace /= xvec_nullspace.norm("l2");

  MatNullSpace petsc_nullspace;
  Vec nvec[1];
  nvec[0] = xvec_nullspace.vec();
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, nvec, &petsc_nullspace);
  MatSetNullSpace(Anest->mat(), petsc_nullspace);

  PETScVector y;
  auto u0tmp1 = std::make_shared<Function>(V);
  auto p0tmp1 = std::make_shared<Function>(Q);
  Anest->init_vectors(y, {u0tmp1->vector(), p0tmp1->vector()});

  PetscBool test_ns = PETSC_TRUE;
  MatNullSpaceTest(petsc_nullspace, Anest->mat(), &test_ns);
  if (test_ns == PETSC_TRUE)
    std::cout << "******* is nullspace: " << test_ns << std::endl;
  else
    std::cout << "******* is NOT nullspace: " << test_ns << std::endl;

  MatSetNullSpace(Anest->mat(), petsc_nullspace);
  // --------------

  // Create near null space basis ---
  VectorSpaceBasis null_space
    = build_nullspace_nested(*V, *u0->vector());

  // Copy vectors
  std::vector<PETScVector> near_nullspace;
  for (std::size_t i = 0; i < null_space.dim(); ++i)
  {
    dolfin_assert(null_space[i]);
    const PETScVector& x
      = null_space[i]->down_cast<PETScVector>();
    near_nullspace.push_back(x);
  }

  // Get pointers to underlying PETSc objects
  std::vector<Vec> petsc_vec(near_nullspace.size());
  for (std::size_t i = 0; i < near_nullspace.size(); ++i)
    petsc_vec[i] = near_nullspace[i].vec();

  MatNullSpace petsc_near_nullspace;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, near_nullspace.size(),
                     petsc_vec.data(), &petsc_near_nullspace);
  MatSetNearNullSpace(A00->mat(), petsc_near_nullspace);

// ---

  // Solve
  PETScKrylovSolver solver;

  // FIXME: Should this be moved into the PETScKrylovSolver  constructor?
  //  KSPSetFromOptions(solver.ksp());
  solver.set_from_options();
  solver.set_operators(Anest, Pnest);

  // Set field-split blocks for preconditioner
  std::vector<la_index> u_dofs, p_dofs;
  Anest->get_block_dofs(u_dofs, 0);
  Anest->get_block_dofs(p_dofs, 1);
  PETScPreconditioner::set_fieldsplit(solver, {u_dofs, p_dofs},
                                              {"0", "1"});
  solver.solve(xvec, bvec);

  // Clean up nullspace
  MatNullSpaceDestroy(&petsc_near_nullspace);
  MatNullSpaceDestroy(&petsc_nullspace);


  XDMFFile xdmf_u(mesh->mpi_comm(), "u0.xdmf");
  xdmf_u.write(*u0);

  XDMFFile xdmf_p(mesh->mpi_comm(), "p0.xdmf");
  xdmf_p.write(*p0);

  double unorm = u0->vector()->norm("l2");
  if (MPI::rank(mesh->mpi_comm()) == 0)
  {
    std::cout << "u norm:" << unorm << std::endl;
  }

  list_timings(TimingClear::clear, {TimingType::wall});

  return 0;
}
