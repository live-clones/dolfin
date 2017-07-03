import dolfin_test.cpp.function
import dolfin_test.function.expression
import dolfin_test.cpp.generation
import dolfin_test.function.functionspace



class MyExpression(dolfin_test.function.expression.UserExpression):
    def eval(self, values, x):
        print("Inside my expression")


mesh = dolfin_test.cpp.generation.UnitSquareMesh(6, 9)
V = dolfin_test.function.functionspace.FunctionSpace(mesh, "Lagrange", 1)
e = MyExpression(V)
print(e.value_rank())

print("---------------------")
mesh = dolfin_test.cpp.generation.UnitCubeMesh(6, 9, 2)
V = dolfin_test.function.functionspace.VectorFunctionSpace(mesh, "Lagrange", 1)
e = MyExpression(V)
print(e.value_rank())
print(e.value_dimension(0))
