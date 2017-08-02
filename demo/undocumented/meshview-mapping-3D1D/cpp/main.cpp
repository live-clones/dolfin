#include <dolfin.h>
#include "MeshView_3D1D.h"

using namespace dolfin;

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    values[0] = 10*exp(-(dx*dx) / 0.02);
  }
};

int main()
{
  // Create mesh and function space
  auto mesh = std::make_shared<UnitCubeMesh>(32, 32, 32);

  EdgeFunction<std::size_t> marker(mesh, 0);
  for (EdgeIterator edge(*mesh); !edge.end(); ++edge)
  {
      auto y = edge->midpoint().y();
      auto z = edge->midpoint().z();
      marker[edge->index()] = (0.5 - DOLFIN_EPS < z && z < 0.5 + DOLFIN_EPS) && (0.5 - DOLFIN_EPS < y && y < 0.5 + DOLFIN_EPS);
  }
  
  std::vector<std::size_t> vertex_map,cell_map;
  auto mapping = std::make_shared<MeshViewMapping>(mesh,vertex_map,cell_map);
  auto submesh = std::make_shared<Mesh>(mapping->create_from_marker(marker, 1));

  // Function spaces associated with each of the function spaces
  auto V1 = std::make_shared<MeshView_3D1D::Form_a00::TestSpace>(mesh); // 3D
  auto V2 = std::make_shared<MeshView_3D1D::Form_a11::TestSpace>(submesh); // 1D

  // Bilinear and linear forms
  MeshView_3D1D::Form_a00 a_3D(V1, V1);
  MeshView_3D1D::Form_a11 a_1D(V2, V2);
  MeshView_3D1D::Form_L0 L_3D(V1);
  MeshView_3D1D::Form_L1 L_1D(V2);

  // Define boundary conditions
  auto zero = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc_3D(V1, zero, boundary);
  DirichletBC bc_1D(V2, zero, boundary);

  // Define RHS
  auto f = std::make_shared<Source>();
  L_3D.f1 = f;
  L_1D.f2 = f;

  // Compute solution
  // Subproblem 3D
  Function u_3D(V1);
  solve(a_3D == L_3D, u_3D, bc_3D);
  // Subproblem 1D
  Function u_1D(V2);
  solve(a_1D == L_1D, u_1D, bc_1D);

  // Save solution in vtk format
  File out_3D("meshview-mapping-3D1D-3Dsol.pvd");
  out_3D << u_3D;
  File out_1D("meshview-mapping-3D1D-1Dsol.pvd");
  out_1D << u_1D;
}


