from ufl import TestFunction, dx

import dolfin_test.cpp.function
import dolfin_test.function.expression
import dolfin_test.cpp.generation
import dolfin_test.function.functionspace

from dolfin_test.cpp.la import EigenVector
from dolfin_test.cpp import MPI
from dolfin_test.cpp.fem import Assembler
from dolfin_test.fem.form import Form
from dolfin_test.function.constant import Constant

class MyExpression(dolfin_test.function.expression.UserExpression):
    def eval(self, values, x):
        print("Inside my expression")


mesh0 = dolfin_test.cpp.generation.UnitSquareMesh(6, 9)
V0 = dolfin_test.function.functionspace.FunctionSpace(mesh0, "Lagrange", 1)
e0 = MyExpression(V0)
print(e0.value_rank())

print("---------------------")
mesh = dolfin_test.cpp.generation.UnitCubeMesh(6, 9, 2)
V = dolfin_test.function.functionspace.VectorFunctionSpace(mesh, "Lagrange", 1)
e = MyExpression(V)
print(e.value_rank())
print(e.value_dimension(0))


e0 = Constant(1.0)
v = TestFunction(V0)
L = e0*v*dx
assembler = Assembler()
b = EigenVector(MPI.comm_world, 0)
assembler.assemble(b, Form(L, [V0]))
