#include <dolfin.h>
#include "MeshView_3D2D.h"

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
    double dy = x[1] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

int main()
{
  // Create mesh and function space
  auto mesh = std::make_shared<UnitCubeMesh>(32, 32, 32);

  FacetFunction<std::size_t> marker(mesh, 0);
  for (FacetIterator facet(*mesh); !facet.end(); ++facet)
  {
      auto x = facet->midpoint().coordinates();
      marker[facet->index()] = 0.5 - DOLFIN_EPS < x[2] && x[2] < 0.5 + DOLFIN_EPS;
  }
  
  std::vector<std::size_t> vertex_map,cell_map;
  auto mapping = std::make_shared<MeshViewMapping>(mesh,vertex_map,cell_map);
  auto submesh = std::make_shared<Mesh>(mapping->create_from_marker(marker, 1));

  // Function spaces associated with each of the function spaces
  auto V1 = std::make_shared<MeshView_3D2D::Form_a00::TestSpace>(mesh); // 3D
  auto V2 = std::make_shared<MeshView_3D2D::Form_a11::TestSpace>(submesh); // 2D

  // Bilinear and linear forms
  MeshView_3D2D::Form_a00 a_3D(V1, V1);
  MeshView_3D2D::Form_a11 a_2D(V2, V2);
  MeshView_3D2D::Form_L0 L_3D(V1);
  MeshView_3D2D::Form_L1 L_2D(V2);

  // Define boundary conditions
  auto zero = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc_3D(V1, zero, boundary);
  DirichletBC bc_2D(V2, zero, boundary);

  // Define RHS
  auto f = std::make_shared<Source>();
  L_3D.f1 = f;
  L_2D.f2 = f;

  // Compute solution
  // Subproblem 3D
  Function u_3D(V1);
  solve(a_3D == L_3D, u_3D, bc_3D);
  // Subproblem 2D
  Function u_2D(V2);
  solve(a_2D == L_2D, u_2D, bc_2D);

  // Save solution in vtk format
  File out_3D("meshview-mapping-3D2D-3Dsol.pvd");
  out_3D << u_3D;
  File out_2D("meshview-mapping-3D2D-2Dsol.pvd");
  out_2D << u_2D;
}


