
import types
from six import string_types
import hashlib

import dolfin_test.cpp as cpp
import dijitso
import ffc

def jit_generate(inside_code, module_name, signature, parameters):

    template_code = """

#include<dolfin.h>
#include<Eigen/Dense>

namespace dolfin
{
  class %(classname)s : public SubDomain
  {
     public:
       %(members)s

       %(classname)s()
          {
            %(constructor)s
          }

       void hello() const
       {
         int a;
         a += 1;
       }

       /// Return true for points inside the sub domain
       bool inside(const Eigen::Ref<Eigen::VectorXd>& x, bool on_boundary) const override
       {
         return %(inside)s;
       }
  };
}

extern "C" __attribute__ ((visibility ("default"))) dolfin::SubDomain * create_%(classname)s()
{
  return new dolfin::%(classname)s;
}

"""
    classname = signature
    code_c = template_code % {"inside": inside_code, "classname": classname, "members": "", "constructor": ""}
    code_h = ""
    depends = []

    print(code_c)

    return code_h, code_c, depends

def compile_subdomain(inside_code):

    module_hash = hashlib.md5(inside_code.encode('utf-8')).hexdigest()
    module_name = "subdomain_" + module_hash
    params = dijitso.params.default_params()
    params['build']['include_dirs'] = ['/usr/include/eigen3',
                                       '/home/chris/src/FEniCS/dolfin/local.garth.feature-pybind11/include',
                                       ffc.backends.ufc.get_include_path()]
    params['build']['libs'] = ['dolfin']
    params['build']['lib_dirs'] = ['/home/chris/src/FEniCS/dolfin/local.garth.feature-pybind11/lib']
    print(params['build'])

    module, signature = dijitso.jit(inside_code, module_name, params,
                                    generate=jit_generate)

    submodule = dijitso.extract_factory_function(module, "create_" + module_name)()
    print("JIT gives:", submodule, module, signature)

    sd = cpp.mesh.make_dolfin_subdomain(submodule)
    return sd

class CompiledSubDomain(cpp.mesh.SubDomain):
    def __init__(self, inside_code):
        self._sd = compile_subdomain(inside_code)
        super().__init__()

    def inside(self, x, on_boundary):
        return self._sd.inside(x, on_boundary)

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

        if not isinstance(args[2], cpp.mesh.SubDomain):
            raise(RuntimeError, "Third argument must be of type SubDomain")
        subdomain = args[2]

        super().__init__(function_space, function, subdomain)
