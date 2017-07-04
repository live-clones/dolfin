# -*- coding: utf-8 -*-
from __future__ import print_function

__all__ = ["Expression"]

# Python imports
import types
from six import add_metaclass
from six import string_types
from six.moves import xrange as range
from functools import reduce
import weakref
import hashlib

import dijitso

# Import UFL and SWIG-generated extension module (DOLFIN C++)
import ufl
from ufl import product
from ufl.utils.indexflattening import flatten_multiindex, shape_to_strides
import dolfin_test.cpp as cpp
import numpy

#from dolfin import warning, error

class UserExpression(ufl.Coefficient, cpp.function.Expression):
    def __init__(self, function_space=None, element=None, degree=None):
        """Create an Expression."""

        #if element is None and degree is None:
        #    raise RuntimeError('UserExpression must specific a FiniteElement or a dgeree')

        #if element is None:
        #    # Create UFL element
        #    element = ufl.FiniteElement(family, mesh.ufl_cell(), degree,
        #                                form_degree=None)

        ufl.Coefficient.__init__(self, function_space)
        print(self.ufl_shape)
        #cpp.function.Expression.__init__(self, 1)
        cpp.function.Expression.__init__(self, self.ufl_shape)

        value_shape = tuple(self.value_dimension(i)
                            for i in range(self.value_rank()))

        #element = _auto_select_element_from_shape(value_shape, degree, cell)
        # Check if scalar, vector or tensor valued
        family = "Lagrange"
        degree = 2
        cell = None
        if len(value_shape) == 0:
            element = ufl.FiniteElement(family, cell, degree)
        elif len(value_shape) == 1:
            element = ufl.VectorElement(family, cell, degree, dim=value_shape[0])
        else:
            element = ufl.TensorElement(family, cell, degree, shape=value_shape)

        # Initialize UFL base class
        ufl_function_space = ufl.FunctionSpace(None, element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())


def jit_generate(statement, module_name, signature, parameters):

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
         values[0] = {statement};
       }}
  }};
}}

extern "C" __attribute__ ((visibility ("default"))) dolfin::Expression * create_{classname}()
{{
  return new dolfin::{classname};
}}

"""

    classname = signature
    code_c = template_code.format(statement=statement, classname=classname,
                                  members= "", constructor="")
    code_h = ""
    depends = []

    return code_h, code_c, depends


def compile_expression(statement):

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

    module_hash = hashlib.md5(statement.encode('utf-8')).hexdigest()
    module_name = "expression_" + module_hash
    module, signature = dijitso.jit(statement, module_name, params,
                                    generate=jit_generate)

    submodule = dijitso.extract_factory_function(module, "create_" + module_name)()
    print("JIT gives:", submodule, module, signature)

    expression = cpp.function.make_dolfin_expression(submodule)
    return expression


class CompiledExpression(ufl.Coefficient, cpp.function.Expression):
    def __new__(cls, statement, degree):
        return compile_expression(statement)

    def __init__(self, statement, degree):
        element = ufl.FiniteElement("Lagrange", None, degree)
        ufl_function_space = ufl.FunctionSpace(None, element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())
