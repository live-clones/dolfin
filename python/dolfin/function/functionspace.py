
import ufl
import dolfin_test.cpp as cpp

class FunctionSpace(ufl.FunctionSpace, cpp.function.FunctionSpace):

    def __init__(self):
        print("hello")
