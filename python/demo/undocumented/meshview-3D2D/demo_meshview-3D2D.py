from dolfin import *

# Create mesh and define function space
mesh = UnitCubeMesh(32, 32, 32)

marker = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
for f in facets(mesh):
    marker[f] = 0.5 - DOLFIN_EPS < f.midpoint().z() < 0.5 + DOLFIN_EPS

submesh = MeshView.create(marker, 1)

V1 = FunctionSpace(mesh, "Lagrange", 1) ## 3D
V2 = FunctionSpace(submesh, "Lagrange", 1) ## 2D

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
# 3D mesh
bc_3D = DirichletBC(V1, u0, boundary)
# 2D submesh
bc_2D = DirichletBC(V2, u0, boundary)

# Define variational problem
# 3D
u_3D = TrialFunction(V1)
v_3D = TestFunction(V1)
# 2D
u_2D = TrialFunction(V2)
v_2D = TestFunction(V2)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2) ) / 0.02)", degree=2)

a_3D = inner(grad(u_3D), grad(v_3D))*dx
L_3D = f*v_3D*dx
u_3D = Function(V1)
solve(a_3D == L_3D, u_3D, bc_3D)

a_2D = inner(grad(u_2D), grad(v_2D))*dx
L_2D = f*v_2D*dx
u_2D = Function(V2)
solve(a_2D == L_2D, u_2D, bc_2D)

# ## Export result
# encoding = XDMFFile.Encoding.HDF5 if has_hdf5() else XDMFFile.Encoding.ASCII

# out_3D = XDMFFile(MPI.comm_world, "meshview-mapping-3D2D-3Dsol.xdmf")
# out_2D = XDMFFile(MPI.comm_world, "meshview-mapping-3D2D-2Dsol.xdmf")

# if MPI.size(MPI.comm_world) > 1 and encoding == XDMFFile.Encoding.ASCII:
#     print("XDMF file output not supported in parallel without HDF5")
# else:
#     out_3D.write(u_3D, encoding)
#     out_2D.write(u_2D, encoding)
