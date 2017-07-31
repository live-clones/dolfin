from dolfin import *
import numpy as np

#class MyNewExpression(dolfin.function.expression.UserExpression):
class MyNewExpression(CompiledExpression):
    #def eval_cell(self, values, x, cell):
    #    print("in eval")
    #    values[0] = 40.0
    def eval(self, values, x):
        print("in my eval")
        values[0] = 20.0

mesh0 = UnitSquareMesh(1, 1)
V0 = FunctionSpace(mesh0, "Lagrange", 1)
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
#e0 = MyNewExpression(V0)
v = TestFunction(V0)
L = e0*v*dx
assembler = Assembler()
b = EigenVector(MPI.comm_world, 0)

#form = Form(L)
#print(type(b), type(form))
assembler.assemble(b, Form(L))

#print(b.array())


#f = CompiledExpression("2.0", degree=1);
#print(dir(f))


#v = np.zeros(1)
#x = np.array([0.0, 0.0])
#f.eval(v, x)
#print(v)
