import dolfin_test.cpp as cpp
from dolfin_test.cpp.generation import UnitSquareMesh
from dolfin_test.function.functionspace import FunctionSpace
from dolfin_test.cpp.fem import Assembler
from dolfin_test.fem.form import Form
from dolfin_test.fem.dirichletbc import DirichletBC, CompiledSubDomain
from dolfin_test.cpp.function import Function, Constant
from dolfin_test.cpp.mesh import SubDomain
from dolfin_test.cpp.la import EigenVector, EigenMatrix, LUSolver
from dolfin_test.cpp import MPI
from dolfin_test.cpp.io import XDMFFile
from ufl import TestFunction, TrialFunction, inner, grad, dx, ds

mesh = UnitSquareMesh(12, 12)
V = FunctionSpace(mesh, "Lagrange", 1)
w = Function(V)

xdmf = XDMFFile("a.xdmf")
xdmf.write(mesh) #, XDMFFile.Encoding.ASCII)
xdmf.write(w) #, XDMFFile.Encoding.ASCII)

DOLFIN_EPS = 1e-14

class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        result = (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS)
        return bool(result)

# boundary = Boundary()

boundary = CompiledSubDomain("x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS")

u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

u = TrialFunction(V)
v = TestFunction(V)
# f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
# g = Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
#L = f*v*dx + g*v*ds
L = v*dx

assembler = Assembler()
A = EigenMatrix()
assembler.assemble(A, Form(a, [V, V]))
print(A.array())

b = EigenVector(MPI.comm_world, 0)
assembler.assemble(b, Form(L, [V]))
bc.apply(b)
bc.apply(A)
# print(b.array())
# print(A.array())

solver = LUSolver(MPI.comm_world, A, "default")

solver.solve(w.vector(), b)

# print (w.vector().array())

# print(w.vector().array())
# solve(a == L, w, bc)

file = XDMFFile("poisson.xdmf")
file.write(w)

# plot(w, interactive=True)
