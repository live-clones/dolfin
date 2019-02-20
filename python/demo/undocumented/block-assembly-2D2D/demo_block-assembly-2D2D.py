from dolfin import *

#Create mesh and define function space
mesh = UnitSquareMesh(30, 30)

marker = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
for c in cells(mesh):
    marker[c] = c.midpoint().y() < 0.5
    #marker[c] = c.midpoint().x() < 0.5

submesh1 = MeshView.create(marker, 1)
submesh2 = MeshView.create(marker, 0)

# Define Dirichlet boundary
def boundarySub1(x):
    return x[1] < DOLFIN_EPS
    #return x[0] < DOLFIN_EPS

def boundarySub2(x):
    return x[1] > 1.0 - DOLFIN_EPS
    #return x[0] > 1.0 - DOLFIN_EPS

#element2D = FiniteElement("Lagrange", triangle, 1)
W1 = FunctionSpace(submesh1, "Lagrange", 1)
W2 = FunctionSpace(submesh2, "Lagrange", 1)
# Define the product function space
V = FunctionSpaceProduct( W1, W2 )

# Define boundary conditions
u0 = Constant(0.0)
# Subdomain 1
bc1 = DirichletBC(V.sub_space(0), u0, boundarySub1)
# Subdomain 2
bc2 = DirichletBC(V.sub_space(1), u0, boundarySub2)
bcs = [bc1,bc2]

# Define variational problem
# Use directly TrialFunction and TestFunction on the product space
# Subdomain 1
u1 = TrialFunction(W1)
v1 = TestFunction(W1)
# Subdomain 2
u2 = TrialFunction(W2)
v2 = TestFunction(W2)
# Mixed
(u1_m,u2_m) = TrialFunction(V)
(v1_m,v2_m) = TestFunction(V)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)

# Weak form - 1
a1 = inner(grad(u1), grad(v1))*dx
L1 = f*v1*dx
u1 = Function(W1)
solve(a1 == L1, u1, bc1)
# Weak form - 1
a2 = inner(grad(u2), grad(v2))*dx
L2 = f*v2*dx
u2 = Function(W2)
solve(a2 == L2, u2, bc2)

# Mixed weak form
a = inner(grad(u1_m), grad(v1_m))*dx + inner(grad(u2_m), grad(v2_m))*dx
L = f*v1_m*dx + f*v2_m*dx
# Solve the problem
sol = Function(V)
solve(a == L, sol, bcs, solver_parameters={"linear_solver":"direct"})

sol1 = sol.sub(0)
sol2 = sol.sub(1)

assert len(u1.vector()) == len(sol1.vector())
for i in range(len(u1.vector())):
    assert abs(sol1.vector()[i] - u1.vector()[i]) < 1e-10
    
assert len(u2.vector()) == len(sol2.vector())
for i in range(len(u2.vector())):
    assert abs(sol2.vector()[i] - u2.vector()[i]) < 1e-10

## Export result
encoding = XDMFFile.Encoding.HDF5 if has_hdf5() else XDMFFile.Encoding.ASCII

out_u1 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-subdomain1-ref.xdmf")
out_u2 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-subdomain2-ref.xdmf")
out_sub1 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-subdomain1.xdmf")
out_sub2 = XDMFFile(MPI.comm_world, "block-assembly-2D2D-subdomain2.xdmf")

if MPI.size(MPI.comm_world) > 1 and encoding == XDMFFile.Encoding.ASCII:
    print("XDMF file output not supported in parallel without HDF5")
else:
    out_u1.write(u1, encoding)
    out_u2.write(u2, encoding)
    out_sub1.write(sol1, encoding)
    out_sub2.write(sol2, encoding)
