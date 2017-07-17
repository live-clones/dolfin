from dolfin import *

# Create mesh
mesh = UnitSquareMesh(6, 9)

# Create function space
V = FunctionSpace(mesh, "Lagrange", 1)

# Variation problem
u = TrialFunction(V)
v = TestFunction(V)

c = Constant(1.0)

a = dot(grad(u), grad(v))*dx
L = c*v*dx

# Assemble
A = assemble(a)
print(A)
