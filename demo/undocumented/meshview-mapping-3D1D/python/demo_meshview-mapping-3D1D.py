from dolfin import *

# Create mesh and define function space
mesh = UnitCubeMesh(32, 32, 32)

marker = EdgeFunction("size_t", mesh, 0)
for e in edges(mesh):
    marker[e] = 0.5 - DOLFIN_EPS < e.midpoint().z() < 0.5 + DOLFIN_EPS and 0.5 - DOLFIN_EPS < e.midpoint().y() < 0.5 + DOLFIN_EPS

submesh = MeshViewMapping.create_from_marker(marker, 1)

V1 = FunctionSpace(mesh, "Lagrange", 1) ## 3D
V2 = FunctionSpace(submesh, "Lagrange", 1) ## 1D

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
# 3D mesh
bc_3D = DirichletBC(V1, u0, boundary)
# 1D submesh
bc_1D = DirichletBC(V2, u0, boundary)

# Define variational problem
# 3D
u_3D = TrialFunction(V1)
v_3D = TestFunction(V1)
# 1D
u_1D = TrialFunction(V2)
v_1D = TestFunction(V2)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) ) / 0.02)", degree=2)

a_3D = inner(grad(u_3D), grad(v_3D))*dx
L_3D = f*v_3D*dx
u_3D = Function(V1)
solve(a_3D == L_3D, u_3D, bc_3D)

a_1D = inner(grad(u_1D), grad(v_1D))*dx
L_1D = f*v_1D*dx
u_1D = Function(V2)
solve(a_1D == L_1D, u_1D, bc_1D)

# Save solution in XDMF format if available
out_3D = XDMFFile(mesh.mpi_comm(), "meshview-mapping-3D1D-3Dsol.xdmf")
out_1D = XDMFFile(mesh.mpi_comm(), "meshview-mapping-3D1D-1Dsol.xdmf")
if has_hdf5():
    out_3D.write(u_3D)
    out_1D.write(u_1D)
else:
    # Save solution in vtk format
    out_3D = File("meshview-mapping-3D1D-3Dsol.pvd")
    out_3D << u_3D
    out_1D = File("meshview-mapping-3D1D-1Dsol.pvd")
    out_1D << u_1D
