from dolfin import *

# Create mesh and define function space
mesh = UnitCubeMesh(32, 32, 32)

marker = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
for f in facets(mesh):
    marker[f] = 0.5 - DOLFIN_EPS < f.midpoint().z() < 0.5 + DOLFIN_EPS

submesh = MeshView.create(marker, 1)

W1 = FunctionSpace(mesh, "Lagrange", 1) ## 3D
W2 = FunctionSpace(submesh, "Lagrange", 1) ## 2D

# Define the product function space
V = FunctionSpaceProduct( W1, W2 )

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
# 3D mesh
bc_3D = DirichletBC(V.sub_space(0), u0, boundary)
# 2D submesh
bc_2D = DirichletBC(V.sub_space(1), u0, boundary)
bcs = [bc_3D,bc_2D]

# Define variational problem
# 3D
u_3D = TrialFunction(W1)
v_3D = TestFunction(W1)
# 2D
u_2D = TrialFunction(W2)
v_2D = TestFunction(W2)
# Use directly TrialFunction and TestFunction on the product space
(u_3D_m,u_2D_m) = TrialFunctions(V)
(v_3D_m,v_2D_m) = TestFunctions(V)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2) ) / 0.02)", degree=2)

# Weak form - 1 (3D problem)
a_3D = inner(grad(u_3D), grad(v_3D))*dx
L_3D = f*v_3D*dx
u_3D = Function(W1)
solve(a_3D == L_3D, u_3D, bc_3D)
# Weak form - 2 (2D problem)
a_2D = inner(grad(u_2D), grad(v_2D))*dx
L_2D = f*v_2D*dx
u_2D = Function(W2)
solve(a_2D == L_2D, u_2D, bc_2D)

# Global weak form (Mixed problem)
a = inner(grad(u_3D_m), grad(v_3D_m))*dx + inner(grad(u_2D_m), grad(v_2D_m))*dx
L = f*v_3D_m*dx + f*v_2D_m*dx

# Solve the problem
sol = Function(V)

# rtol=1e-4
# print("******************************************************************")
# print("NOTE : Relative tolerance for Krylov solver has been set to", rtol)
# print("With default solver/preconditioner, convergence is very slow.")
# print("We need to find a better configuration.")
# print("******************************************************************")
# solve(a == L, sol, bcs=bcs, solver_parameters={"krylov_solver":{"relative_tolerance":rtol, "maximum_iterations":10000}})

# Direct solve
solve(a == L, sol, bcs, solver_parameters={"linear_solver":"direct"})
sol_3D = sol.sub(0)
sol_2D = sol.sub(1)

assert len(u_3D.vector()) == len(sol_3D.vector())
for i in range(len(u_3D.vector())):
    assert abs(u_3D.vector()[i] - sol_3D.vector()[i]) < 1e-10
    
assert len(u_2D.vector()) == len(sol_2D.vector())
for i in range(len(u_2D.vector())):
    assert abs(u_2D.vector()[i] - sol_2D.vector()[i]) < 1e-10

## Export result
encoding = XDMFFile.Encoding.HDF5 if has_hdf5() else XDMFFile.Encoding.ASCII

out_u3D = XDMFFile(MPI.comm_world, "block-assembly-3D2D-3Dsol-ref.xdmf")
out_u2D = XDMFFile(MPI.comm_world, "block-assembly-3D2D-2Dsol-ref.xdmf")
out_3D = XDMFFile(MPI.comm_world, "block-assembly-3D2D-3Dsol.xdmf")
out_2D = XDMFFile(MPI.comm_world, "block-assembly-3D2D-2Dsol.xdmf")

if MPI.size(MPI.comm_world) > 1 and encoding == XDMFFile.Encoding.ASCII:
    print("XDMF file output not supported in parallel without HDF5")
else:
    out_u3D.write(u_3D, encoding)
    out_u2D.write(u_2D, encoding)
    out_3D.write(sol_3D, encoding)
    out_2D.write(sol_2D, encoding)
