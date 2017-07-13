from dolfin import *

#Create mesh and define function space
mesh = UnitSquareMesh(32, 32)

marker = CellFunction("size_t", mesh, 0)
for c in cells(mesh):
    marker[c] = c.midpoint().x() < 0.5

submesh1 = MeshViewMapping.create_from_marker(marker, 1)
submesh2 = MeshViewMapping.create_from_marker(marker, 0)

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

# Define new measures so that we can sort the integrals
a_prod = inner(grad(u1), grad(v1))*dx + inner(grad(u2), grad(v2))*dx
L_prod = f*v1*dx + f*v2*dx
# Subdomain 1
a1 = extract_blocks(a_prod,0,0)
L1 = extract_blocks(L_prod,0)
sol1 = Function(V.sub_space(0))
solve(a1 == L1, sol1, bc1)
# Subdomain 2
a2 = extract_blocks(a_prod,1,1)
L2 = extract_blocks(L_prod,1)
sol2 = Function(V.sub_space(1))
solve(a2 == L2, sol2, bc2)

# Save solution in XDMF format if available
out_sub1 = XDMFFile(mesh.mpi_comm(), "formsplitter-product-subdomain1.xdmf")
out_sub2 = XDMFFile(mesh.mpi_comm(), "formsplitter-product-subdomain2.xdmf")
if has_hdf5():
    out_sub1.write(sol1)
    out_sub2.write(sol2)
else:
    # Save solution in vtk format
    out_sub1 = File("formsplitter-product-subdomain1.pvd")
    out_sub1 << sol1
    out_sub2 = File("formsplittere-product-subdomain2.pvd")
    out_sub2 << sol2
