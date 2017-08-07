# -*- coding: utf-8 -*-

import dolfin.cpp as cpp
from dolfin.function.expression import BaseExpression, _select_element

from six import string_types
import hashlib
import dijitso

def jit_generate(class_data, module_name, signature, parameters):

    template_code = """

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>
#include <dolfin/mesh/MeshFunction.h>
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

  }};
}}

PYBIND11_MODULE({classname}, m)
{{
  py::class_<dolfin::{classname}, std::shared_ptr<dolfin::{classname}>, dolfin::Expression>
    (m, "{classname}", py::dynamic_attr())
    .def(py::init<>())
    {properties}
    .def("eval", &dolfin::{classname}::eval);
}}

"""

    _property_code = """    .def_property("{a}", [](dolfin::{classname}& self){{return self.{a};}},
                       [](dolfin::{classname}& self, {type} val){{self.{a} = val;}})
"""

    statements = class_data["statements"]
    statement_code = ""

    property_values = class_data["properties"]
    for key in property_values:
        statement_code += "         const auto& ref = *{a};\n".format(a=key)

    for i, val in enumerate(statements):
        statement_code += "          values[" + str(i) + "] = " + val + ";\n"


    classname = signature
    members_code = ""
    constructor_code = ""
    properties_code = ""

    for key in property_values:
        if isinstance(property_values[key], float):
            members_code += "double {a};\n".format(a=key)
            properties_code += _property_code.format(a=key, type='double', classname=classname)
        else:
            # *** temporary hacked in MeshFunction
            members_code += "std::shared_ptr<dolfin::MeshFunction<std::size_t>> {a};\n".format(a=key)
            properties_code += _property_code.format(a=key, type='std::shared_ptr<dolfin::MeshFunction<std::size_t>>', classname=classname)


    # Set the value_shape
    if len(statements) > 1:
        constructor_code += "_value_shape.push_back(" + str(len(statements)) + ");"

    code_c = template_code.format(statement=statement_code, properties=properties_code,
                                  members=members_code, constructor=constructor_code,
                                  classname=classname)
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
    # FIXME: need to locate Python libs and pybind11 - this works on my Mac...
    from distutils import sysconfig
    import pybind11
    params = dijitso.params.default_params()
    pyversion = "python" + sysconfig.get_config_var("LDVERSION")
    params['cache']['lib_prefix'] = ""
    params['cache']['lib_basename'] = ""
    params['cache']['lib_loader'] = "import"
    params['build']['include_dirs'] = d["include_dirs"] + [pybind11.get_include(), sysconfig.get_config_var("INCLUDEDIR") + "/" + pyversion]
    params['build']['libs'] = d["libraries"] + [ pyversion ]
    params['build']['lib_dirs'] = d["library_dirs"] + [sysconfig.get_config_var("LIBDIR")]
    params['build']['cxxflags'] += ('-fno-lto',)

    print(params)

    if isinstance(statements, string_types):
        statements = tuple((statements,))

    if not isinstance(statements, tuple):
        raise RuntimeError("Expression must be a string, or a tuple of strings")

    keystr = tuple(key for key in properties)
    #FIXME: include properties (class variable names) in hash
    module_hash = hashlib.md5("".join(statements + keystr).encode('utf-8')).hexdigest()
    module_name = "dolfin_expression_" + module_hash

    class_data = {"statements": statements, "properties": properties}

    module, signature = dijitso.jit(class_data,
                                    module_name, params,
                                    generate=jit_generate)

    print('pybind11 JIT: ', module)

    expression = getattr(module, signature)()

    # Set values
    for key in properties:
        setattr(expression, key, properties[key])

    return expression

class CompiledExpressionPyBind11(BaseExpression):
    """JIT Expressions"""

    def __init__(self, cpp_code=None, *args, **kwargs):

        # Remove arguments that are used in Expression creation
        element = kwargs.pop("element", None)
        degree = kwargs.pop("degree", None)
        cell = kwargs.pop("cell", None)
        domain = kwargs.pop("domain", None)
        name = kwargs.pop("name", None)
        label = kwargs.pop("label", None)
        mpi_comm = kwargs.pop("mpi_comm", None)

        # Save properties for checking later
        self._properties = kwargs
        for k in self._properties:
            if not isinstance(k, string_types):
                raise KeyError("Invalid key")

        if cpp_code is not None:
            self._cpp_object = compile_expression(cpp_code, self._properties)

        # Deduce element type if not provided
        if element is None:
            if degree is None:
                raise KeyError("Must supply element or degree")
            value_shape = tuple(self.value_dimension(i)
                                for i in range(self.value_rank()))
            element = _select_element(family=None, cell=None, degree=degree,
                                      value_shape=value_shape)

        BaseExpression.__init__(self, cell=cell, element=element, domain=domain,
                                name=name, label=label)

    # This is added dynamically in the intialiser to allow checking of
    # eval in user classes.
    def __getattr__(self, name):
        "Pass attributes through to (JIT compiled) Expression object"
        if name in self._properties.keys():
            return getattr(self._cpp_object, name)
        else:
            raise(AttributeError)

    def __setattr__(self, name, value):
        if name.startswith("_"):
            super().__setattr__(name, value)
        elif name in self._properties.keys():
            setattr(self._cpp_object, name, value)


