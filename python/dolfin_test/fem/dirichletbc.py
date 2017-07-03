
import types
from six import string_types

import dolfin_test.cpp as cpp
import dijitso


class CompiledSubdomain(cpp.mesh.SubDomain):
    def inside(self, x, on_boundary):
        pass


def jit_generate(compile_string, module_name, signature, parameters):

    code_c = compile_string
    code_h = ""
    depends = []

    return code_h, code_c, depends

def compiled_subdomain(inside_code):

    classname="subdomain"

    template_code = """
#include<Eigen/Dense>

namespace dolfin
{
  class %(classname)s: public SubDomain
  {
     public:
       %(members)s

       %(classname)s()
          {
            %(constructor)s
          }

       /// Return true for points inside the sub domain
       bool inside(const Eigen::VectorXd& x, bool on_boundary) const
       {
         return %(inside)s;
       }
  };
}
""" % {"classname": classname, "members": "", "constructor": "",
       "inside": inside_code}

    module_name = "subdomain"
    params = None

    module, signature = dijitso.jit(template_code, module_name, params,
                                    generate=jit_generate)

    print(module, signature)

    return None


class DirichletBC(cpp.fem.DirichletBC):
    def __init__(self, *args, **kwargs):

        if len(args) != 3:
            raise(RuntimeError, "Not yet supported")

        if not isinstance(args[0], cpp.function.FunctionSpace):
            raise(RuntimeError)
        function_space = args[0]

        if not isinstance(args[1], cpp.function.GenericFunction):
            raise(RuntimeError)
        function = args[1]

        if isinstance(args[2], cpp.mesh.SubDomain):
            subdomain = args[2]
        elif isinstance(args[2], string_types):
            subdomain = compiled_subdomain(args[2])
        else:
            raise(RuntimeError)

        cpp.fem.DirichletBC.__init__(self, function_space, function, subdomain)
