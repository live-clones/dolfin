from dolfin import *

# Sub domains for Dirichlet boundary conditions
class DirichletBoundary1(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 0.0) < DOLFIN_EPS and on_boundary
class DirichletBoundary2(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary

def newton_solver_parameters():
    return{"nonlinear_solver": "newton",
           "newton_solver": {"linear_solver": "gmres"}}

# Create meshes
n = 16
mesh = UnitSquareMesh(n, n)

marker = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
for c in cells(mesh):
    marker[c] = c.midpoint().x() < 0.5

submesh1 = MeshView.create(marker, 1)
submesh2 = MeshView.create(marker, 0)

# Create function spaces
W1 = FunctionSpace(submesh1, "Lagrange", 1)
W2 = FunctionSpace(submesh2, "Lagrange", 1)
# Define the mixed function space W = W1 x W2
W = MixedFunctionSpace( W1, W2 )

# Define boundary conditions
g = Constant(1.0)
bc1 = DirichletBC(W1, g, DirichletBoundary1())
bc2 = DirichletBC(W2, g, DirichletBoundary2())
bcs = [bc1,bc2]

f = Expression("x[0]*sin(x[1])", degree=2)

# Define mixed-domains variational problem
(v1,v2) = TestFunctions(W)
u = Function(W)
u1 = u.sub(0)
u2 = u.sub(1)

dx1 = Measure("dx", domain=W.sub_space(0).mesh())
dx2 = Measure("dx", domain=W.sub_space(1).mesh()) 

F1 = inner((1 + u1**2)*grad(u1), grad(v1))*dx1 - f*v1*dx1
F2 = inner((1 + u2**2)*grad(u2), grad(v2))*dx2 - f*v2*dx2
F = F1 + F2

# Compute solution - ref monodomain problems
solve(F1 == 0, u1, bc1)
u1_ref = u1.copy(deepcopy=True)
solve(F2 == 0, u2, bc2)
u2_ref = u2.copy(deepcopy=True)

# Compute solution - mixed-domains problem
solve(F == 0, u, bcs, solver_parameters=newton_solver_parameters())
# solve(F == 0, u, bcs, solver_parameters={"nonlinear_solver":"snes"}) # Not available yet

for i in range(len(u1.vector())):
    assert abs(u1_ref.vector()[i] - u1.vector()[i]) < 1e-8
for i in range(len(u2.vector())):
    assert abs(u2_ref.vector()[i] - u2.vector()[i]) < 1e-8

## Export result
encoding = XDMFFile.Encoding.HDF5 if has_hdf5() else XDMFFile.Encoding.ASCII
out_ref1 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-nonlinear-subdomain1-ref.xdmf")
out_ref2 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-nonlinear-subdomain2-ref.xdmf")
out_sub1 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-nonlinear-subdomain1.xdmf")
out_sub2 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-nonlinear-subdomain2.xdmf")

if MPI.size(MPI.comm_world) > 1 and encoding == XDMFFile.Encoding.ASCII:
    print("XDMF file output not supported in parallel without HDF5")
else:
    out_ref1.write(u1_ref, encoding)
    out_ref2.write(u2_ref, encoding)
    out_sub1.write(u1, encoding)
    out_sub2.write(u2, encoding)
