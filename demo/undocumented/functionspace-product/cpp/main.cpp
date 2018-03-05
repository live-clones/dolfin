#include <dolfin.h>
#include "FunctionSpaceProduct.h"

using namespace dolfin;

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

class DirichletBoundarySubdomain1 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS;
  }
};

class DirichletBoundarySubdomain2 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] > 1.0 - DOLFIN_EPS;
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
  // Create mesh
  auto mesh = std::make_shared<UnitSquareMesh>(32, 32);

  // Try to create two mesh views from this 2D mesh
  MeshFunction<std::size_t> marker(mesh, mesh->topology().dim(), 0);
  for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      auto x = cell->midpoint().coordinates();
      marker[cell->index()] = x[0] < 0.5;
    }

  std::vector<std::size_t> vertex_map,cell_map;
  auto mapping = std::make_shared<MeshView>(mesh,vertex_map,cell_map);
  auto submesh1 = std::make_shared<Mesh>( mapping->create(marker,1) );
  auto submesh2 = std::make_shared<Mesh>( mapping->create(marker,0) );

  // Function spaces associated with each of the function spaces
  auto V1=std::make_shared<FunctionSpaceProduct::Form_a00::TestSpace>( submesh1 );
  auto V2=std::make_shared<FunctionSpaceProduct::Form_a11::TestSpace>( submesh2 );

  // Bilinear and linear forms (defined in the ufl file)
  FunctionSpaceProduct::Form_a00 a1(V1, V1);
  FunctionSpaceProduct::Form_a11 a2(V2, V2);
  FunctionSpaceProduct::Form_L0 L1(V1);
  FunctionSpaceProduct::Form_L1 L2(V2);

  // Define boundary conditions
  auto zero = std::make_shared<Constant>(0.0);
  auto boundarySubdomain1=std::make_shared<DirichletBoundarySubdomain1>();
  auto boundarySubdomain2=std::make_shared<DirichletBoundarySubdomain2>();
  DirichletBC bc1(V1, zero, boundarySubdomain1);
  DirichletBC bc2(V2, zero, boundarySubdomain2);

  // Define RHS
  //auto f = std::make_shared<Constant>(0.0);
  auto f = std::make_shared<Source>();
  L1.f1=f;
  L2.f2=f;

  // Compute solution
  // Subproblem 1
  Function u1(V1);
  solve(a1 == L1, u1, bc1);
  // Subproblem 1
  Function u2(V2);
  solve(a2 == L2, u2, bc2);

  // Save solution in vtk format
  File out_sub1("functionspace-product-subdomain1.pvd");
  out_sub1 << u1;
  File out_sub2("functionspace-product-subdomain2.pvd");
  out_sub2 << u2;
}