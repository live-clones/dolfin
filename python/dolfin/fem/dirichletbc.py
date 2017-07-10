
import types
from six import string_types

import dolfin.cpp as cpp
from dolfin.mesh.subdomain import CompiledSubDomain

class DirichletBC(cpp.fem.DirichletBC):
    def __init__(self, *args, **kwargs):

        if len(args) != 3:
            raise(RuntimeError, "Not yet supported")

        if not isinstance(args[0], cpp.function.FunctionSpace):
            raise(RuntimeError, "First argument must be of type FunctionSpace")
        function_space = args[0]

        if not isinstance(args[1], cpp.function.GenericFunction):
            raise(RuntimeError, "Second argument must be of type GenericFunction")
        function = args[1]

        if isinstance(args[2], cpp.mesh.SubDomain):
            subdomain = args[2]
        elif isinstance(args[2], string_types):
            subdomain = CompiledSubDomain(args[2])
        else:
            raise(RuntimeError, "Third argument must be of type SubDomain or string")

        super().__init__(function_space, function, subdomain)
