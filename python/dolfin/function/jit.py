# -*- coding: utf-8 -*-
import hashlib
from six import string_types
import dijitso
import dolfin.cpp as cpp


def jit_generate(class_data, module_name, signature, parameters):
    """TODO: document"""

    template_code = """

#include <dolfin/function/Expression.h>
#include <Eigen/Dense>

namespace dolfin
{{
  class {classname} : public Expression
  {{
     public:
       {members}

       {classname}()
          {{
            {constructor}
          }}

       void eval(Eigen::Ref<Eigen::VectorXd> values, const Eigen::Ref<Eigen::VectorXd> x) const override
       {{
{statement}
       }}

       void set_property(std::string name, double value)
       {{
{set_props}
       }}

       double get_property(std::string name) const
       {{
{get_props}
       }}

  }};
}}

extern "C" __attribute__ ((visibility ("default"))) dolfin::Expression * create_{classname}()
{{
  return new dolfin::{classname};
}}

"""
    _get_props = """          if (name == "{name}") return {name};"""
    _set_props = """          if (name == "{name}") {{ {name} = value; return; }}"""

    statements = class_data["statements"]
    statement = ""
    for i, val in enumerate(statements):
        statement += "          values[" + str(i) + "] = " + val + ";\n"

    # Set the value_shape
    constructor = ""
    if len(statements) > 1:
        constructor += "_value_shape.push_back(" + str(len(statements)) + ");"

    members = ""
    set_props = ""
    get_props = ""

    # Add code for setting and getting property values
    properties = class_data["properties"]
    for k in properties:
        value = properties[k]
        members += "double " + k + ";\n"
        set_props += _set_props.format(name=k)
        get_props += _get_props.format(name=k)

    classname = signature
    code_c = template_code.format(statement=statement, classname=classname,
                                  members=members, constructor=constructor,
                                  set_props=set_props, get_props=get_props)
    code_h = ""
    depends = []

    return code_h, code_c, depends


def compile_expression(statements, properties):
    """Compile a user C(++) string to a Python object"""

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

    if not isinstance(statements, (string_types, tuple, list)):
        raise RuntimeError("Expression must be a string, or a list or tuple of strings")

    if isinstance(statements, string_types):
        statements = (statements,)

    class_data = {'statements': statements, 'properties': properties}

    hash_str = str(statements)
    module_hash = hashlib.md5(hash_str.encode('utf-8')).hexdigest()
    module_name = "dolfin_expression_" + module_hash

    try:
        module, signature = dijitso.jit(class_data, module_name, params,
                                        generate=jit_generate)
        submodule = dijitso.extract_factory_function(module, "create_" + module_name)()
    except:
        raise RuntimeError("Unable to compile C++ code with dijitso")

    expression = cpp.function.make_dolfin_expression(submodule)

    # Set properties to initial values
    for k in properties:
        expression.set_property(k, properties[k])

    return expression
