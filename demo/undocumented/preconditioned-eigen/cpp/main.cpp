// Copyright (C) 2014 Garth N. Wells
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

#include <dolfin.h>
#include "Elasticity.h"

using namespace dolfin;

class Left : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS;
  }
};


int main()
{
  // Create Mesh and
  UnitCubeMesh mesh(1, 1, 1);
  Elasticity::FunctionSpace V(mesh);

  Left left;
  Constant zero(0.0, 0.0, 0.0);
  DirichletBC bc(V, zero, left);

  // Set elasticity parameters
  double E  = 10.0;
  double nu = 0.3;
  Constant mu(E / (2*(1 + nu)));
  Constant lambda(E*nu / ((1 + nu)*(1 - 2*nu)));

  // Define variational problem
  Elasticity::BilinearForm a(V, V);
  a.mu = mu; a.lmbda = lambda;

  Elasticity::LinearForm L(V);
  Constant f(1.0, 0.0, 0.0);
  L.f = f;

  PETScMatrix A;
  PETScVector b;
  SystemAssembler assembler(a, L, bc);
  assembler.assemble(A, b);

  PETScVector x(b);
  x.zero();

  SLEPcEigenSolver eigen_solver(A);
  eigen_solver.parameters["solver"] = "lapack";
  std::cout << "Solving eigenproblem" << std::endl;
  eigen_solver.solve();
  std::cout << "Done solving eigenproblem" << std::endl;
  double lr, lc;
  for (std::size_t i = 0; i < A.size(0); ++i)
  {
    eigen_solver.get_eigenvalue(lr, lc, i);
    std::cout << lr << ", " << lc  << std::endl;
  }

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCILU);

  KSPSetOperators(ksp, A.mat(), A.mat());
  KSPSetFromOptions(ksp);

  /*
  KSPSolve(ksp, b.vec(), x.vec());
  PetscInt size = 0;
  std::vector<double> r(A.size(0)), c(A.size(0));
  KSPComputeEigenvaluesExplicitly(ksp, r.size(), r.data(), c.data());
  std::sort(r.begin(), r.end(), std::greater<double>());
  std::cout << "Eigenvalues (explicit computation): " << std::endl;
  for (std::size_t i = 0; i < r.size(); ++i)
    std::cout << " " << r[i] << ", " << std::endl;
  */

  EPS eps;
  EPSCreate(PetscObjectComm((PetscObject)ksp), &eps);
  Mat Amat, Pmat;
  KSPGetOperators(ksp,  &Amat, &Pmat);
  EPSSetOperators(eps, Amat, Pmat);

  ST st;
  KSP stksp;
  EPSGetST(eps, &st);
  STGetKSP(st, &stksp);
  if (stksp)
  {
    std::cout << "***Set pc " << std::endl;
    PC pc;
    KSPGetPC(ksp, &pc);
    KSPSetPC(stksp, pc);
    KSPSetType(stksp, KSPPREONLY);
  }

  EPSSetDimensions(eps, A.size(0), PETSC_DECIDE, PETSC_DECIDE);

  EPSSetFromOptions(eps);

  EPSSolve(eps);

  int nconv = 0;
  EPSGetConverged(eps, &nconv);
  std::cout << "SLEPC, num converged eigenvalues: " << nconv << "  (" << A.size(0) << ")" << std::endl;
  for (int j = 0; j < nconv; j++)
  {
    double r, c;
    Vec xr, xc;
    EPSGetEigenpair(eps, j, &r, &c, xr, xc);
    std::cout << r << ", " << c << std::endl;
  }

 return 0;
}
