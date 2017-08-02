import types
from six import string_types
import hashlib

import dolfin.cpp as cpp
import dijitso
import ffc


def jit_generate(class_data, module_name, signature, parameters):

    template_code = """
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/common/Array.h>
#include <dolfin/math/basic.h>
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

       // Return true for points inside the sub domain
       bool inside(const Eigen::Ref<Eigen::VectorXd> x, bool on_boundary) const final
       {{
         return {inside};
       }}

       void set_property(std::string name, double value)
       {{
{set_props}
       }}

       double get_property(std::string name) const
       {{
{get_props}
         return 0.0;
       }}

  }};
}}

extern "C" DLL_EXPORT dolfin::SubDomain * create_{classname}()
{{
  return new dolfin::{classname};
}}

"""
    _set_prop = """ if (name == "{name}") {name} = value;\n"""
    _get_prop = """ if (name == "{name}") return {name};\n"""

    print('Class data = ', class_data)
    inside_code = class_data['inside_code']

    members = ""
    get_props = ""
    set_props = ""
    for k in class_data['properties']:
        members += " double " + k + ";\n"
        get_props += _get_prop.format(name = k)
        set_props += _set_prop.format(name = k)

    classname = signature
    code_c = template_code.format(inside=inside_code, classname=classname,
                                  members=members, constructor="", get_props=get_props, set_props=set_props)
    code_h = ""
    depends = []

    return code_h, code_c, depends


def compile_subdomain(inside_code, properties):

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

    class_data = {'inside_code': inside_code, 'properties': properties}

    hash_str = inside_code + str(properties.keys())
    module_hash = hashlib.md5(hash_str.encode('utf-8')).hexdigest()
    module_name = "dolfin_subdomain_" + module_hash

    try:
        module, signature = dijitso.jit(class_data, module_name, params,
                                        generate=jit_generate)
        submodule = dijitso.extract_factory_function(module, "create_" + module_name)()
    except:
        raise RuntimeError("Unable to compile C++ code with dijitso")

    sub_domain = cpp.mesh.make_dolfin_subdomain(submodule)

    for k in properties:
        sub_domain.set_property(k, properties[k])

    return sub_domain

class CompiledSubDomain(cpp.mesh.SubDomain):
    def __new__(cls, inside_code, **kwargs):
        properties = kwargs
        return compile_subdomain(inside_code, properties)
