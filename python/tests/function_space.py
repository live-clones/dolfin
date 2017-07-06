import ufl
import dolfin.function.functionspace
import dolfin.function.constant
import dolfin.cpp.generation

# Create mesh
mesh = dolfin.cpp.generation.UnitSquareMesh(6, 9)

# Create function space
V = dolfin.function.functionspace.FunctionSpace(mesh, "Lagrange", 1)

# Variation problem
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

c = dolfin.function.constant.Constant(1.0)

a = ufl.dot(ufl.grad(u), ufl.grad(v))*(ufl.dx)
L = c*v*(ufl.dx)

# Assemble
