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
# Use directly TrialFunction and TestFunction on the product space
(u_3D,u_2D) = TrialFunction(V)
(v_3D,v_2D) = TestFunction(V)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2) ) / 0.02)", degree=2)

# Global weak form
a = inner(grad(u_3D), grad(v_3D))*dx + inner(grad(u_2D), grad(v_2D))*dx
L = f*v_3D*dx + f*v_2D*dx

# Solve the problem
sol = Function(V)
rtol=1e-4
print("******************************************************************")
print("NOTE : Relative tolerance for Krylov solver has been set to", rtol)
print("With default solver/preconditioner, convergence is very slow.")
print("We need to find a better configuration.")
print("******************************************************************")
solve(a == L, sol, bcs=bcs, solver_parameters={"krylov_solver":{"relative_tolerance":rtol, "maximum_iterations":10000}})

# extract components of the solution
sol_3D = sol.sub(0)
sol_2D = sol.sub(1)

## Save solution in vtk format
out_3D = File("block-assembly-3D2D-3Dsol.pvd")
out_3D << sol_3D
out_2D = File("block-assembly-3D2D-2Dsol.pvd")
out_2D << sol_2D
