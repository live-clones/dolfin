//
// Demo for MatNest and DOLFIN
// solving Stokes on a 2D square
//
// Chris Richardson Jan 2016

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

// Tests that the constant pressure field is in the nullspace
void test_pressure_nullspace(std::shared_ptr<dolfin::Mesh> mesh)
{
  class Boundary : public dolfin::SubDomain
  {
    bool inside(const dolfin::Array<double>& x, bool on_boundary) const
    { return on_boundary; }
  };

  // Create function spaces
  auto V = std::make_shared<Stokes::Form_a00::TestSpace>(mesh);
  auto Q = std::make_shared<Stokes::Form_a11::TestSpace>(mesh);

  // Velocity BC
  auto boundary = std::make_shared<Boundary>();
  auto zero_vector = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  auto bc = std::make_shared<dolfin::DirichletBC>(V, zero_vector, boundary);

  // Create forms
  auto a00 = std::make_shared<Stokes::Form_a00>(V, V);
  auto a01 = std::make_shared<Stokes::Form_a01>(Q, V);
  auto a10 = std::make_shared<Stokes::Form_a10>(V, Q);

  // Junk we don't need. Need to fix DOLFIN interface
  auto L0 = std::make_shared<Stokes::Form_L0>(V);
  auto L1 = std::make_shared<Stokes::Form_L1>(Q);
  L0->f = std::make_shared<dolfin::Constant>(0.0, 0.0, 0.0);
  L1->g = std::make_shared<dolfin::Constant>(0.0);

  auto b0 = std::make_shared<dolfin::PETScVector>();
  auto b1 = std::make_shared<dolfin::PETScVector>();

  // Assemble all blocks
  dolfin::SystemAssembler assemblerA0({a00, a01, a10, NULL}, {L0, L1},
                                      {{bc}, {}});
  auto A00 = std::make_shared<dolfin::PETScMatrix>();
  auto A01 = std::make_shared<dolfin::PETScMatrix>();
  auto A10 = std::make_shared<dolfin::PETScMatrix>();
  assemblerA0.assemble({A00, A01, A10, NULL}, {b0, b1});

  // Create pressure nullspace vector
  auto u = std::make_shared<dolfin::Function>(V);
  auto p = std::make_shared<dolfin::Function>(Q);

  dolfin::Constant one(1.0);
  dolfin::Constant one_vec(0.0, 0.0, 0.0);
  p->interpolate(one);
  u->interpolate(one_vec);
  double norm = p->vector()->norm("l2");
  (*p->vector()) /= norm;

  // Create nested matrix
  dolfin::PETScNestMatrix A({A00, A01, A10, NULL});

  // Create nullspace vector
  dolfin::PETScVector nullspace;
  A.init_vectors(nullspace, {u->vector(), p->vector()});

  // Create PETSc nullspace
  MatNullSpace petsc_nullspace;
  PetscErrorCode ierr;
  Vec nvec[1];
  nvec[0] = nullspace.vec();
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1,
                            nvec, &petsc_nullspace);
  MatSetNullSpace(A.mat(), petsc_nullspace);

  // Test nullspace
  PetscBool test_ns = PETSC_TRUE;
  MatNullSpaceTest(petsc_nullspace, A.mat(), &test_ns);
  if (test_ns == PETSC_TRUE)
    std::cout << "*** Is nullspace: " << test_ns << std::endl;
  else
    std::cout << "**Is NOT nullspace: " << test_ns << std::endl;

  MatNullSpaceDestroy(&petsc_nullspace);
}

// Function to compute the near nullspace A00 operator
dolfin::VectorSpaceBasis build_nullspace_nested(const dolfin::FunctionSpace& V,
                                                const dolfin::GenericVector& x)
{
  // Get subspaces
  auto V0 = V.sub(0);
  auto V1 = V.sub(1);
  auto V2 = V.sub(2);

  // Create vectors for nullspace basis
  std::vector<std::shared_ptr<dolfin::GenericVector>> basis(3);
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
  dolfin::VectorSpaceBasis vector_space(basis);
  vector_space.orthonormalize();
  return vector_space;
}


std::array<std::shared_ptr<dolfin::Function>, 2>
  solve_nested(std::shared_ptr<dolfin::Mesh> mesh,
               std::shared_ptr<dolfin::SubDomain> top_bottom,
               std::shared_ptr<dolfin::SubDomain> left_edge,
               std::shared_ptr<dolfin::GenericFunction> flow_velocity,
               std::shared_ptr<dolfin::GenericFunction> pleft,
               std::shared_ptr<dolfin::GenericFunction> f,
               std::shared_ptr<dolfin::GenericFunction> zero,
               bool pc_operator)
{
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
  auto bcV = std::make_shared<dolfin::DirichletBC>(V, flow_velocity,
                                                   top_bottom);

  // Pressure BC
  //auto bcQ = std::make_shared<dolfin::DirichletBC>(Q, pleft, left_edge);

  // Assemble all blocks with BCs
  dolfin::SystemAssembler assemblerA0({a00, a01, a10, a11}, {L0, L1},
                                      {{bcV}, {}});
  auto A00 = std::make_shared<dolfin::PETScMatrix>();
  auto A01 = std::make_shared<dolfin::PETScMatrix>();
  auto b0 = std::make_shared<dolfin::PETScVector>();
  auto A10 = std::make_shared<dolfin::PETScMatrix>();
  auto b1 = std::make_shared<dolfin::PETScVector>();

  auto A11 = std::make_shared<dolfin::PETScMatrix>();
  std::cout << "Assemble A00, A01, A10, A11\n";
  assemblerA0.assemble({A00, A01, A10, A11}, {b0, b1});

  // Assemble a mass matrix (QxQ) for the preconditioner
  auto one = std::make_shared<dolfin::Constant>(1.0);
  a11->c = one;
  dolfin::SystemAssembler assemblerP(a11, L1, {});
  auto P11 = std::make_shared<dolfin::PETScMatrix>();
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
  std::vector<std::shared_ptr<const dolfin::GenericMatrix>> Amats, Pmats;
  std::shared_ptr<dolfin::PETScNestMatrix> Anest, Pnest;
  if (pc_operator)
  {
    Amats = {A00, A01, A10, NULL};
    Pmats = {A00, NULL, NULL, P11};
    Pnest = std::make_shared<dolfin::PETScNestMatrix>(Pmats);
  }
  else
  {
    Amats = {A00, A01, A10, A11};
  }
  Anest = std::make_shared<dolfin::PETScNestMatrix>(Amats);

  // Create VecNest vectors for RHS and solution vector
  auto u0 = std::make_shared<dolfin::Function>(V);
  auto p0 = std::make_shared<dolfin::Function>(Q);
  dolfin::PETScVector xvec, bvec;
  Anest->init_vectors(xvec, {u0->vector(), p0->vector()});
  Anest->init_vectors(bvec, {b0, b1});

  // ---------
  // Create pressure nullspace vector
  auto u0tmp = std::make_shared<dolfin::Function>(V);
  auto p0tmp = std::make_shared<dolfin::Function>(Q);
  p0tmp->interpolate(*one);

  dolfin::Constant one_vec(0.0, 0.0, 0.0);
  u0tmp->interpolate(one_vec);

  dolfin::PETScVector xvec_nullspace;
  Anest->init_vectors(xvec_nullspace, {u0tmp->vector(), p0tmp->vector()});
  xvec_nullspace /= xvec_nullspace.norm("l2");

  MatNullSpace petsc_nullspace;
  Vec nvec[1];
  nvec[0] = xvec_nullspace.vec();
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, nvec, &petsc_nullspace);
  MatSetNullSpace(Anest->mat(), petsc_nullspace);

  dolfin::PETScVector y;
  auto u0tmp1 = std::make_shared<dolfin::Function>(V);
  auto p0tmp1 = std::make_shared<dolfin::Function>(Q);
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
  dolfin::VectorSpaceBasis null_space
    = build_nullspace_nested(*V, *u0->vector());

  // Copy vectors
  std::vector<dolfin::PETScVector> near_nullspace;
  for (std::size_t i = 0; i < null_space.dim(); ++i)
  {
    dolfin_assert(null_space[i]);
    const dolfin::PETScVector& x
      = null_space[i]->down_cast<dolfin::PETScVector>();
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
  dolfin::PETScKrylovSolver solver;

  // FIXME: Should this be moved into the PETScKrylovSolver  constructor?
  KSPSetFromOptions(solver.ksp());

  if (pc_operator)
    solver.set_operators(Anest, Pnest);
  else
    solver.set_operator(Anest);

  // Set field-split blocks for preconditioner
  std::vector<dolfin::la_index> u_dofs, p_dofs;
  Anest->get_block_dofs(u_dofs, 0);
  Anest->get_block_dofs(p_dofs, 1);
  dolfin::PETScPreconditioner::set_fieldsplit(solver, {u_dofs, p_dofs},
                                              {"0", "1"});
  solver.solve(xvec, bvec);

  // Test solver on A00 block.
  /*
  {
    dolfin::PETScKrylovSolver solver_ref;
    dolfin::PETScOptions::set("pc_type", "gamg");
    solver_ref.set_operator(A00);
    solver_ref.solve(*(u0->vector()), *b0);
  }
  */

  // Clean up nullspace
  MatNullSpaceDestroy(&petsc_near_nullspace);
  MatNullSpaceDestroy(&petsc_nullspace);

  return {u0, p0};
}


// Function for no-slip boundary condition for velocity
class Source : public dolfin::Expression
{
public:

  Source() : Expression(3) {}

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0] = x[1];
    values[1] = x[1]*x[2];
    values[2] = 0.0;
  }
};

// Sub domain for top and bottom
class TopBottom : public dolfin::SubDomain
{
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const
  {
    return on_boundary;
    //return std::abs(1.0 - x[1]) < DOLFIN_EPS || std::abs(x[1]) < DOLFIN_EPS;
  }
};

class LeftEdge : public dolfin::SubDomain
{
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS;
  }
};


using namespace dolfin;

int main(int argc, char *argv[])
{
  // SCOTCH is crashing on ARCHER
  parameters["mesh_partitioner"] = "ParMETIS";
  parameters["timer_prefix"] = "TICKS:";

  // Application parameters
  Parameters application_parameters("application_parameters");
  application_parameters.add("size", 12);
  application_parameters.add("mode", 0);

  // Parse application parameters
  application_parameters.parse(argc, argv);
  const int size = application_parameters["size"];
  const int mode = application_parameters["mode"];

  std::cout << "mode = " << mode << "\n";

// Parse PETSc parameters
  parameters.parse(argc, argv);

  //PETScOptions::set("ksp_view");
  PETScOptions::set("ksp_monitor_true_residual");

  // Create mesh
  auto mesh = std::make_shared<UnitCubeMesh>(size, size, size);

  // Velocity BC
  auto top_bottom = std::make_shared<TopBottom>();
  auto flow_velocity = std::make_shared<Constant>(0.0, 0.0, 0.0) ;

  // Pressure BC
  auto left_edge = std::make_shared<LeftEdge>();
  auto pleft = std::make_shared<Constant>(12.0);

  //auto f = std::make_shared<dolfin::Constant>(0.5, 1.0, 1.0);
  auto f = std::make_shared<Source>();
  auto zero = std::make_shared<dolfin::Constant>(0.0);

  // Common fieldsplit options
  dolfin::PETScOptions::set("ksp_type", "minres");
  //dolfin::PETScOptions::set("ksp_type", "gmres");
  dolfin::PETScOptions::set("pc_type", "fieldsplit");

  bool pc_operator = true;
  {
    dolfin::PETScOptions::set("ksp_view");

    dolfin::PETScOptions::set("pc_fieldsplit_type", "additive");

    // AMG
    dolfin::PETScOptions::set("fieldsplit_0_ksp_type", "preonly");
    //dolfin::PETScOptions::set("fieldsplit_0_pc_type", "lu");
    //dolfin::PETScOptions::set("fieldsplit_0_pc_factor_mat_solver_package",
    //                          "mumps");

    // "gamg" or "hypre"
    std::string mg = "hypre";
    if (mg == "gamg")
    {
      dolfin::PETScOptions::set("fieldsplit_0_pc_type", "gamg");
      dolfin::PETScOptions::set("fieldsplit_0_pc_gamg_coarse_eq_limit", 10000);

      dolfin::PETScOptions::set("fieldsplit_0_pc_gamg_threshold", 0.05);
      //dolfin::PETScOptions::set("fieldsplit_0_pc_gamg_square_graph", 2);
      dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_type", "chebyshev");
      dolfin::PETScOptions::set("fieldsplit_0_mg_levels_pc_type", "sor");
      dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_max_it", 3);

      dolfin::PETScOptions::set("fieldsplit_0_mg_levels_esteig_ksp_type", "cg");
      dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_chebyshev_esteig_steps", 20);
      dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_chebyshev_esteig_random");

      dolfin::PETScOptions::set("fieldsplit_0_mg_coarse_ksp_type", "preonly");
      dolfin::PETScOptions::set("fieldsplit_0_mg_coarse_pc_type", "lu");
      dolfin::PETScOptions::set("fieldsplit_0_mg_coarse_pc_factor_mat_solver_package", "mumps");
    }
    else if (mg == "hypre")
    {
      dolfin::PETScOptions::set("fieldsplit_0_pc_type", "hypre");
      dolfin::PETScOptions::set("fieldsplit_0_hypre_type", "boomeramg");
      dolfin::PETScOptions::set("fieldsplit_1_pc_hypre_boomeramg_strong_threshold", 0.5);
    }
    else
      dolfin::error("Unknown A00 pc");

    dolfin::PETScOptions::set("fieldsplit_1_ksp_type", "preonly");
    dolfin::PETScOptions::set("fieldsplit_1_pc_type", "hypre");
    dolfin::PETScOptions::set("fieldsplit_1_hypre_type", "boomeramg");
    dolfin::PETScOptions::set("fieldsplit_1_pc_hypre_boomeramg_strong_threshold", 0.5);

//dolfin::PETScOptions::set("fieldsplit_1_pc_type", "lu");
    //dolfin::PETScOptions::set("fieldsplit_1_pc_factor_mat_solver_package",
    //                          "mumps");

    // LU
    /*
    dolfin::PETScOptions::set("fieldsplit_0_ksp_type", "preonly");
    dolfin::PETScOptions::set("fieldsplit_0_pc_type", "lu");
    dolfin::PETScOptions::set("fieldsplit_0_pc_factor_mat_solver_package",
                              "mumps");

    dolfin::PETScOptions::set("fieldsplit_1_ksp_type", "preonly");
    dolfin::PETScOptions::set("fieldsplit_1_pc_type", "lu");
    dolfin::PETScOptions::set("fieldsplit_1_pc_factor_mat_solver_package",
                              "mumps");
    */
  }

  // Exact Schur complement preconditioner. Should converge in three
  // iterations.
  /*
  bool pc_operator = false;
  dolfin::PETScOptions::set("pc_type", "fieldsplit");
  dolfin::PETScOptions::set("pc_fieldsplit_type", "schur");
  dolfin::PETScOptions::set("pc_fieldsplit_schur_fact_type", "diag");

  dolfin::PETScOptions::set("fieldsplit_0_ksp_type", "preonly");
  dolfin::PETScOptions::set("fieldsplit_0_ksp_type", "cg");
  //dolfin::PETScOptions::set("fieldsplit_0_ksp_rtol", 1.0e-9);
  dolfin::PETScOptions::set("fieldsplit_0_pc_type", "gamg");
  dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_type", "chebyshev");
  dolfin::PETScOptions::set("fieldsplit_0_mg_levels_pc_type", "jacobi");
  dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_max_it", 2);
  dolfin::PETScOptions::set("fieldsplit_0_mg_levels_esteig_ksp_type", "cg");
  dolfin::PETScOptions::set("fieldsplit_0_mg_levels_ksp_chebyshev_esteig_steps",
                            50);

  //dolfin::PETScOptions::set("fieldsplit_0_ksp_type", "preonly");
  //dolfin::PETScOptions::set("fieldsplit_0_pc_type", "lu");
  //dolfin::PETScOptions::set("fieldsplit_0_pc_factor_mat_solver_package",
  //                          "superlu_dist");

  // Use tight tolerance on Schur complement solver to get 'near'
  // exact solve
  dolfin::PETScOptions::set("fieldsplit_1_ksp_type", "bcgs");
  dolfin::PETScOptions::set("fieldsplit_1_ksp_rtol", 1.0e-12);
  dolfin::PETScOptions::set("fieldsplit_1_pc_type", "lsc");
  //dolfin::PETScOptions::set("fieldsplit_1_lsc_pc_type", "ilu");
  //dolfin::PETScOptions::set("fieldsplit_1_ksp_monitor_short");
  dolfin::PETScOptions::set("fieldsplit_1_ksp_constant_null_space");
  dolfin::PETScOptions::set("fieldsplit_1_lsc_ksp_constant_null_space");
  */
  //test_pressure_nullspace(mesh);

  double unorm0(-1.0), unorm1(-1.0), unorm2(-1.0);

  // Solve using PETSc MatNest
  if (MPI::rank(mesh->mpi_comm()) == 0)
    std::cout << "Solving MatNest\n";
  auto soln0 = solve_nested(mesh, top_bottom, left_edge, flow_velocity,
                            pleft, f, zero, pc_operator);

  unorm0 = soln0[0]->vector()->norm("l2");
  dolfin::XDMFFile xdmf_u(mesh->mpi_comm(), "u0.xdmf");
  xdmf_u.write(*soln0[0]);
  //dolfin::XDMFFile xdmf_p(mesh->mpi_comm(), "p0.xdmf");
  //xdmf_p.write(*soln0[1]);

  if (dolfin::MPI::rank(mesh->mpi_comm()) == 0)
  {
    std::cout << "u norm (mat nest):         " << unorm0 << std::endl;
  }

  dolfin::list_timings(dolfin::TimingClear::clear, {dolfin::TimingType::wall});

  return 0;
}
