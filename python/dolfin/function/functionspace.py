
import ufl
import dolfin_test.cpp as cpp

class FunctionSpace(ufl.FunctionSpace, cpp.function.FunctionSpace):

    def __init__(self, *args, **kwargs):
        print('Intialise FunctionSpace with: ', args)

        print(len(args))

        if len(args) == 3:

            mesh, family, degree = args
            element = ufl.FiniteElement(family, mesh.cell_type(), degree, form_degree=None)

            print(element)
