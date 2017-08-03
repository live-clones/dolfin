
import types
from six import string_types

import ufl
import dolfin.cpp as cpp
from dolfin.mesh.subdomain import CompiledSubDomain
from dolfin.function.constant import Constant
from dolfin.fem.projection import project

class AutoSubDomain(cpp.mesh.SubDomain):
    "Wrapper class for creating a SubDomain from an inside() function."

    def __init__(self, inside_function):
        "Create SubDomain subclass for given inside() function"

        # Check that we get a function
        if not isinstance(inside_function, types.FunctionType):
            raise RuntimeError("bcs.py",
                               "auto-create subdomain",
                               "Expecting a function (not %s)" %
                               str(type(inside_function)))
        self.inside_function = inside_function

        # Check the number of arguments
        if inside_function.__code__.co_argcount not in (1, 2):
            raise RuntimeError("bcs.py",
                               "auto-create subdomain",
                               "Expecting a function of the form inside(x) or inside(x, on_boundary)")
        self.num_args = inside_function.__code__.co_argcount

        cpp.mesh.SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        "Return true for points inside the subdomain"

        if self.num_args == 1:
            return self.inside_function(x)
        else:
            return self.inside_function(x, on_boundary)


class DirichletBC(cpp.fem.DirichletBC):
    def __init__(self, *args, **kwargs):

        method = kwargs.pop("method", None)

        if len(args) != 3:
            raise RuntimeError("Not yet supported")

        # Special case for value specified as float, tuple or similar
        if len(args) >= 2 and not isinstance(args[1], cpp.function.GenericFunction):
            if isinstance(args[1], ufl.classes.Expr):
                expr = project(args[1], args[0])
            else:
                expr = Constant(args[1])  # let Constant handle all problems
            args = args[:1] + (expr,) + args[2:]

        # Special case for sub domain specified as a function
        if len(args) >= 3 and isinstance(args[2], types.FunctionType):
            sub_domain = AutoSubDomain(args[2])
            args = args[:2] + (sub_domain,) + args[3:]

        if not isinstance(args[0], cpp.function.FunctionSpace):
            raise RuntimeError("First argument must be of type FunctionSpace")
        V = args[0]

        if isinstance(args[1], float) or isinstance(args[1], int):
            u = cpp.function.Constant(float(args[1]))
        elif isinstance(args[1], ufl.Coefficient):
            # Extract cpp object
            u = args[1].cpp_object()
        elif isinstance(args[1], cpp.function.GenericFunction):
            u = args[1]
        elif not isinstance(args[1], cpp.function.GenericFunction):
            raise RuntimeError("Second argument must be of type GenericFunction: ", args[1], type(args[1]))

        if isinstance(args[2], cpp.mesh.SubDomain):
             subdomain = args[2]
        elif isinstance(args[2], string_types):
             subdomain = CompiledSubDomain(args[2])
        else:
            raise RuntimeError("Third argument must be of type SubDomain or string")

        super().__init__(V, u, subdomain)
