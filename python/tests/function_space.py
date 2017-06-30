import ufl
import dolfin_test.function.functionspace
import dolfin_test.cpp.generation

# Create mesh
mesh = dolfin_test.cpp.generation.UnitSquareMesh(6, 9)

# Create function space
V = dolfin_test.function.functionspace.FunctionSpace(mesh, "Lagrange", 1)

# Variation problem
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

a = ufl.dot(ufl.grad(u), ufl.grad(v))*(ufl.dx)
L = 1.0*v*(ufl.dx)

# Assemble

# TODO
