#include <dolfin.h>
#include "MeshView_3D3D.h"

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
    double dz = x[2] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy + dz*dz) / 0.02);
  }
};

int main()
{
  // Create mesh and function space
  auto mesh = std::make_shared<UnitCubeMesh>(32, 32, 32);
  auto V = std::make_shared<MeshView_3D3D::FunctionSpace>(mesh);

  CellFunction<std::size_t> marker(mesh, 0);
  for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      auto x = cell->midpoint().coordinates();
      marker[cell->index()] = x[0] < 0.5;
    }

  std::vector<std::size_t> vertex_map,cell_map;
  auto mapping = std::make_shared<MeshViewMapping>(mesh,vertex_map,cell_map);
  auto submesh1 = std::make_shared<Mesh>( mapping->create_from_marker(marker,1) );
  auto submesh2 = std::make_shared<Mesh>( mapping->create_from_marker(marker,0) );

  // Function spaces associated with each of the function spaces
  auto V1=std::make_shared<MeshView_3D3D::FunctionSpace>( submesh1 );
  auto V2=std::make_shared<MeshView_3D3D::FunctionSpace>( submesh2 );

  // Bilinear and linear forms
  MeshView_3D3D::BilinearForm a(V, V);
  MeshView_3D3D::LinearForm L(V);

  MeshView_3D3D::BilinearForm a1(V1, V1);
  MeshView_3D3D::BilinearForm a2(V2, V2);
  MeshView_3D3D::LinearForm L1(V1);
  MeshView_3D3D::LinearForm L2(V2);

  // Define boundary conditions
  auto zero = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, zero, boundary);

  auto boundarySubdomain1=std::make_shared<DirichletBoundarySubdomain1>();
  auto boundarySubdomain2=std::make_shared<DirichletBoundarySubdomain2>();
  DirichletBC bc1(V1, zero, boundarySubdomain1);
  DirichletBC bc2(V2, zero, boundarySubdomain2);

  // Define RHS
  auto f = std::make_shared<Source>();
  L.f=f;

  L1.f=f;
  L2.f=f;

  // Compute solution
  // Global problem
  Function u(V);
  solve(a == L, u, bc);
  // Subproblem 1
  Function u1(V1);
  solve(a1 == L1, u1, bc1);
  // Subproblem 1
  Function u2(V2);
  solve(a2 == L2, u2, bc2);

  // Save solution in XDMF format if available
  XDMFFile out_global(mesh->mpi_comm(), "meshview-mapping-3D3D-global.xdmf");
  XDMFFile out_sub1(mesh->mpi_comm(), "meshview-mapping-3D3D-subdomain1.xdmf");
  XDMFFile out_sub2(mesh->mpi_comm(), "meshview-mapping-3D3D-subdomain2.xdmf");
  if( has_hdf5() )
    {
      out_global.write(u);
      out_sub1.write(u1);
      out_sub2.write(u2);
    }
  else
    {
      // Save solution in vtk format
      File out_global("meshview-mapping-3D3D-global.pvd");
      out_global << u;
      File out_sub1("meshview-mapping-3D3D-subdomain1.pvd");
      out_sub1 << u1;
      File out_sub2("meshview-mapping-3D3D-subdomain2.pvd");
      out_sub2 << u2;
    }
}


