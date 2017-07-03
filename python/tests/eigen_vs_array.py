
import timeit

setup="""
from dolfin import UnitSquareMesh, FunctionSpace, Function, SubDomain, Constant, DirichletBC, DOLFIN_EPS

mesh = UnitSquareMesh(100,100)
V = FunctionSpace(mesh, "CG", 3)
f = Function(V)

class Boundary(SubDomain):
    def inside(self, x, inside):
        return (x[0] < DOLFIN_EPS)

boundary = Boundary()
u0 = Constant(1.0)
bc = DirichletBC(V, u0, boundary)
"""

t=timeit.timeit('bc.apply(f.vector())', setup=setup, number=400)
print(t)

setup="""
from dolfin_test.cpp.generation import UnitSquareMesh
from dolfin_test.function.functionspace import FunctionSpace
from dolfin_test.fem.dirichletbc import DirichletBC
from dolfin_test.cpp.function import Function, Constant
from dolfin_test.cpp.mesh import SubDomain


mesh = UnitSquareMesh(100,100)
V = FunctionSpace(mesh, "CG", 3)
f = Function(V)

DOLFIN_EPS=1e-14
class Boundary(SubDomain):
    def inside(self, x, inside):
        return bool(x[0] < DOLFIN_EPS)

boundary = Boundary()
u0 = Constant(1.0)
bc = DirichletBC(V, u0, boundary)
"""

t=timeit.timeit('bc.apply(f.vector())', setup=setup, number=400)
print(t)
