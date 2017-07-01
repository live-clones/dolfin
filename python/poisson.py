
from dolfin_test.cpp.generation import UnitSquareMesh
from dolfin_test.function.functionspace import FunctionSpace
from dolfin_test.cpp.function import Function, Constant
from dolfin_test.cpp.fem import DirichletBC
from dolfin_test.cpp.mesh import SubDomain

mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)
w = Function(V)
print(w.vector())

DOLFIN_EPS = 1e-14

class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

boundary = Boundary()
print(boundary)

u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)
quit()

u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds


solve(a == L, w, bc)

file = File("poisson.pvd")
file << w

# plot(w, interactive=True)
