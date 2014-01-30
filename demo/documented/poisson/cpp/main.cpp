// Copyright (C) 2006-2011 Anders Logg
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
// First added:  2006-02-07
// Last changed: 2013-03-11
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 10*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)
//
// and boundary conditions given by
//
//     u(x, y) = 0        for x = 0 or x = 1
// du/dn(x, y) = sin(5*x) for y = 0 or y = 1

#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(5*x[0]);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};



int main(int argc, char* argv[])
{
  SlepcInitialize(&argc,&argv, (char*)0, (char*)0);


  //Parameters parameters = GlobalParameters();
  //parameters.parse(argc, argv);

  // Create mesh and function space
  UnitSquareMesh mesh(6, 6);
  Poisson::FunctionSpace V(mesh);

  // Define boundary condition
  Constant u0(0.0);
  DirichletBoundary boundary;
  DirichletBC bc(V, u0, boundary);

  // Define variational forms
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);

  Source f;
  dUdN g;
  L.f = f;
  L.g = g;

  // Compute solution
  Function u(V);
  //solve(a == L, u, bc);

  PETScMatrix _A;
  PETScVector _b;
  SystemAssembler assembler(a, L, bc);
  assembler.assemble(_A, _b);

  PETScVector _x(_b);
  _x.zero();

  Mat A = _A.mat();
  Vec b = _b.vec();
  Vec x = _x.vec();


  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A, SAME_PRECONDITIONER);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, b, x);

  PetscInt size = 0;
  VecGetSize(x, &size);
  std::vector<double> r(size), c(size);
  KSPComputeEigenvaluesExplicitly(ksp, r.size(), r.data(), c.data());
  std::cout << "Eigenvalues (explicit computation): " << std::endl;
  for (std::size_t i = 0; i < r.size(); ++i)
    std::cout << " " << r[i] << ", " << c[i] << std::endl;

  EPS eps;
  EPSCreate(PetscObjectComm((PetscObject)ksp), &eps);
  Mat Amat, Pmat;
  KSPGetOperators(ksp, &Amat, &Pmat, NULL);
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

  //EPSSetDimensions(eps, size/4, PETSC_DECIDE, PETSC_DECIDE);

  EPSSetFromOptions(eps);

  EPSSolve(eps);

  int nconv = 0;
  EPSGetConverged(eps, &nconv);
  std::cout << "SLEPC, num converged eigenvalues: " << nconv << std::endl;
  for (int j = 0; j < nconv; j++)
  {
    double r, c;
    Vec xr, xc;
    EPSGetEigenpair(eps, j, &r, &c, xr, xc);
    std::cout << r << ", " << c << std::endl;
  }

  std::cout << "done" << std::endl;

  // Save solution in VTK format
  //File file("poisson.pvd");
  //f/ile << u;

  // Plot solution
  //plot(u);
  //interactive();

  return 0;
}
