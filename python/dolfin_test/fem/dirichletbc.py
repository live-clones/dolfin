import types
from six import string_types
import hashlib

import dolfin_test.cpp as cpp
import dijitso
import ffc


def jit_generate(inside_code, module_name, signature, parameters):

    template_code = """

#include <dolfin/mesh/SubDomain.h>
#include <Eigen/Dense>

namespace dolfin
{{
  class {classname} : public SubDomain
  {{
     public:
       {members}

       {classname}()
          {{
            {constructor}
          }}

       void hello() const
       {{
         int a;
         a += 1;
       }}

       /// Return true for points inside the sub domain
       bool inside(const Eigen::Ref<Eigen::VectorXd>& x, bool on_boundary) const override
       {{
         return {inside};
       }}
  }};
}}

extern "C" __attribute__ ((visibility ("default"))) dolfin::SubDomain * create_{classname}()
{{
  return new dolfin::{classname};
}}

"""

    classname = signature
    code_c = template_code.format(inside=inside_code, classname=classname,
                                  members= "", constructor="")
    print(code_c)
    code_h = ""
    depends = []

    return code_h, code_c, depends


def compile_subdomain(inside_code):

    import pkgconfig
    if not pkgconfig.exists('dolfin'):
        raise RuntimeError("Could not find DOLFIN pkg-config file. Please make sure appropriate paths are set.")

    # Get pkg-config data
    d = pkgconfig.parse('dolfin')

    # Set compiler/build options
    params = dijitso.params.default_params()
    params['build']['include_dirs'] = d["include_dirs"]
    params['build']['libs'] = d["libraries"]
    params['build']['lib_dirs'] = d["library_dirs"]

    module_hash = hashlib.md5(inside_code.encode('utf-8')).hexdigest()
    module_name = "subdomain_" + module_hash
    module, signature = dijitso.jit(inside_code, module_name, params,
                                    generate=jit_generate)

    submodule = dijitso.extract_factory_function(module, "create_" + module_name)()
    print("JIT gives:", submodule, module, signature)

    sub_domain = cpp.mesh.make_dolfin_subdomain(submodule)
    return sub_domain


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
