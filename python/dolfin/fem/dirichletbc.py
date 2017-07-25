
import types
from six import string_types

import ufl
import dolfin.cpp as cpp
from dolfin.mesh.subdomain import CompiledSubDomain

class DirichletBC(cpp.fem.DirichletBC):
    def __init__(self, *args, **kwargs):

        method = kwargs.pop("method", None)

        if len(args) != 3:
            raise RuntimeError("Not yet supported")

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
            raise RuntimeError("Second argument must be of type GenericFunction")

        if isinstance(args[2], cpp.mesh.SubDomain):
             subdomain = args[2]
        elif isinstance(args[2], string_types):
             subdomain = CompiledSubDomain(args[2])
        else:
            raise RuntimeError("Third argument must be of type SubDomain or string")

        super().__init__(V, u, subdomain)
