import ufl
import dolfin_test.function.functionspace
import dolfin_test.cpp.generation

mesh = dolfin_test.cpp.generation.UnitSquareMesh(6, 9)

#V = dolfin_test.cpp.function.FunctionSpace()
V = dolfin_test.function.functionspace.FunctionSpace(mesh, "Lagrange", 1)
