# -*- coding: utf-8 -*-
import hashlib
from six import string_types
import dijitso
import dolfin.cpp as cpp

from dolfin.jit.jit import compile_class, _math_header

def jit_generate(class_data, module_name, signature, parameters):
    """TODO: document"""

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

#include <dolfin/function/Expression.h>
#include <Eigen/Dense>

{math_header}

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

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {{
{statement}
       }}

       void set_property(std::string name, double value) override
       {{
{set_props}
       }}

       double get_property(std::string name) const override
       {{
{get_props}
       }}

  }};
}}

extern "C" DLL_EXPORT dolfin::Expression * create_{classname}()
{{
  return new dolfin::{classname};
}}

"""
    _get_props = """          if (name == "{name}") return {name};"""
    _set_props = """          if (name == "{name}") {{ {name} = value; return; }}"""

    statements = class_data["statements"]
    statement = ""
    if isinstance(statements, string_types):
        statement += "          values[0] = " + statements + ";\n"
    else:
        for i, val in enumerate(statements):
            statement += "          values[" + str(i) + "] = " + val + ";\n"

    constructor = ""
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

    # Set the value_shape
    for dim in class_data['value_shape']:
        constructor += "_value_shape.push_back(" + str(dim) + ");"

    classname = signature
    code_c = template_code.format(statement=statement, classname=classname,
                                  members=members, constructor=constructor,
                                  set_props=set_props, get_props=get_props,
                                  math_header=_math_header)
    code_h = ""
    depends = []

    return code_h, code_c, depends


def compile_expression(statements, properties):

    cpp_data = {'statements': statements,
                'properties': properties,
                'name': 'expression',
                'jit_generate': jit_generate}

    expression = compile_class(cpp_data)
    return expression
