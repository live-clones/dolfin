import ufl
import dolfin_test.function.functionspace
import dolfin_test.function.constant
import dolfin_test.cpp.generation

# Create mesh
mesh = dolfin_test.cpp.generation.UnitSquareMesh(6, 9)

# Create function space
V = dolfin_test.function.functionspace.FunctionSpace(mesh, "Lagrange", 1)

# Variation problem
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

c = dolfin_test.function.constant.Constant(1.0)

a = ufl.dot(ufl.grad(u), ufl.grad(v))*(ufl.dx)
L = c*v*(ufl.dx)

# Assemble
