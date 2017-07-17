import numpy as np

from ufl import TestFunction, dx

import dolfin.cpp.function
import dolfin.function.expression
import dolfin.cpp.generation
import dolfin.function.functionspace

from dolfin.cpp.la import EigenVector
from dolfin.cpp import MPI
from dolfin.cpp.fem import Assembler
from dolfin.fem.form import Form
from dolfin.function.constant import Constant
from dolfin.function.expression import CompiledExpression

#class MyNewExpression(dolfin.function.expression.UserExpression):
class MyNewExpression(dolfin.function.expression.UserExpression):
    #def eval_cell(self, values, x, cell):
    #    print("in eval")
    #    values[0] = 40.0
    def eval(self, values, x):
        print("in my eval")
        values[0] = 20.0

mesh0 = dolfin.cpp.generation.UnitSquareMesh(1, 1)
V0 = dolfin.function.functionspace.FunctionSpace(mesh0, "Lagrange", 1)
e0 = MyNewExpression(V0)
print(e0.value_rank())

values = np.zeros(1)
e0.eval(values, (1.0, 1.0))
print(values)
print("---------------------")
#mesh = dolfin.cpp.generation.UnitCubeMesh(6, 9, 2)
#V = dolfin.function.functionspace.VectorFunctionSpace(mesh, "Lagrange", 1)
#e = MyNewExpression(V)
#print(e.value_rank())
#print(e.value_dimension(0))

#e0 = Constant(1.0)
e0 = MyNewExpression(V0)
v = TestFunction(V0)
L = e0*v*dx
assembler = Assembler()
b = EigenVector(MPI.comm_world, 0)

#form = Form(L)
#print(type(b), type(form))
assembler.assemble(b, form)

#print(b.array())


#f = CompiledExpression("2.0", degree=1);
#print(dir(f))


#v = np.zeros(1)
#x = np.array([0.0, 0.0])
#f.eval(v, x)
#print(v)
