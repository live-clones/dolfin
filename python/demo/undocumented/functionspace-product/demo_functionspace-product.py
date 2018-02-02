from dolfin import *

#Create mesh and define function space
mesh = UnitSquareMesh(32, 32)

marker = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
for c in cells(mesh):
    marker[c] = c.midpoint().x() < 0.5

submesh1 = MeshView.create(marker, 1)
submesh2 = MeshView.create(marker, 0)

# Define Dirichlet boundary
def boundarySub1(x):
    return x[0] < DOLFIN_EPS

def boundarySub2(x):
    return x[0] > 1.0 - DOLFIN_EPS

#element2D = FiniteElement("Lagrange", triangle, 1)
W1 = FunctionSpace(submesh1, "Lagrange", 1)
W2 = FunctionSpace(submesh2, "Lagrange", 2)

# Define the product function space
V = FunctionSpaceProduct( W1, W2 )

# Define boundary conditions
u0 = Constant(0.0)
# Subdomain 1
bc1 = DirichletBC(V.sub_space(0), u0, boundarySub1)
# Subdomain 2
bc2 = DirichletBC(V.sub_space(1), u0, boundarySub2)

# Define variational problem
# Use directly TrialFunction and TestFunction on the product space
(u1,u2) = TrialFunction(V)
(v1,v2) = TestFunction(V)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)

# Subdomain 1
a1 = inner(grad(u1), grad(v1))*dx
L1 = f*v1*dx
sol1 = Function(V.sub_space(0))
solve(a1 == L1, sol1, bc1)
# Subdomain 2
a2 = inner(grad(u2), grad(v2))*dx
L2 = f*v2*dx
sol2 = Function(V.sub_space(1))
solve(a2 == L2, sol2, bc2)

# Save solution in vtk format
out_sub1 = File("functionspace-product-subdomain1.pvd")
out_sub1 << sol1
out_sub2 = File("functionspace-product-subdomain2.pvd")
out_sub2 << sol2
