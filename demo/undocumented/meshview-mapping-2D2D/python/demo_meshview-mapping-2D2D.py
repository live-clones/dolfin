from dolfin import *

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

marker = CellFunction("size_t", mesh, 0)
for c in cells(mesh):
    marker[c] = c.midpoint().x() < 0.5

submesh1 = MeshViewMapping.create_from_marker(marker, 1)
submesh2 = MeshViewMapping.create_from_marker(marker, 0)

V1 = FunctionSpace(submesh1, "Lagrange", 1)
V2 = FunctionSpace(submesh2, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

def boundarySub1(x):
    return x[0] < DOLFIN_EPS

def boundarySub2(x):
    return x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
# Whole domain
bc = DirichletBC(V, u0, boundary)
# Subdomain 1
bc1 = DirichletBC(V1, u0, boundarySub1)
# Subdomain 2
bc2 = DirichletBC(V2, u0, boundarySub2)

# Define variational problem
# Whole domain
u = TrialFunction(V)
v = TestFunction(V)
# Subdomain 1
u1 = TrialFunction(V1)
v1 = TestFunction(V1)
# Subdomain 2
u2 = TrialFunction(V2)
v2 = TestFunction(V2)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)

# Whole domain
a = inner(grad(u), grad(v))*dx
L = f*v*dx
u = Function(V)
solve(a == L, u, bc)
# Subdomain 1
a1 = inner(grad(u1), grad(v1))*dx
L1 = f*v1*dx
u1 = Function(V1)
solve(a1 == L1, u1, bc1)
# Subdomain 2
a2 = inner(grad(u2), grad(v2))*dx
L2 = f*v2*dx
u2 = Function(V2)
solve(a2 == L2, u2, bc2)

# Save solution in XDMF format if available
out_global = XDMFFile(mesh.mpi_comm(), "meshview-mapping-2D2D-global.xdmf")
out_sub1 = XDMFFile(mesh.mpi_comm(), "meshview-mapping-2D2D-subdomain1.xdmf")
out_sub2 = XDMFFile(mesh.mpi_comm(), "meshview-mapping-2D2D-subdomain2.xdmf")
if has_hdf5():
    out_global.write(u)
    out_sub1.write(u1)
    out_sub2.write(u2)
else:
    # Save solution in vtk format
    out_global = File("meshview-mapping-2D2D-global.pvd")
    out_global << u
    out_sub1 = File("meshview-mapping-2D2D-subdomain1.pvd")
    out_sub1 << u1
    out_sub2 = File("meshview-mapping-2D2D-subdomain2.pvd")
    out_sub2 << u2

