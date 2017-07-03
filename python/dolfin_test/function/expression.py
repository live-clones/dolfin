# -*- coding: utf-8 -*-
"""This module handles the Expression class in Python.

The Expression class needs special handling and is not mapped directly
by SWIG from the C++ interface. Instead, a new Expression class is
created which inherits both from the DOLFIN C++ Expression class and
the ufl Coefficient class.

The resulting Expression class may thus act both as a variable in a
UFL form expression and as a DOLFIN C++ Expression.

This module make heavy use of creation of Expression classes and
instantiation of these dynamically at runtime.

The whole logic behind this somewhat magic behaviour is handle by:

  1) function __new__ in the Expression class
  2) meta class ExpressionMetaClass
  3) function compile_expressions from the compiledmodule/expression
     module
  4) functions create_compiled_expression_class and
     create_python_derived_expression_class

The __new__ method in the Expression class take cares of the logic
when the class Expression is used to create an instance of Expression,
see use cases 1-4 in the docstring of Expression.

The meta class ExpressionMetaClass take care of the logic when a user
subclasses Expression to create a user-defined Expression, see use
cases 3 in the docstring of Expression.

The function compile_expression is a JIT compiler. It compiles and
returns different kinds of cpp.Expression classes, depending on the
arguments. These classes are sent to the
create_compiled_expression_class.

"""

# Copyright (C) 2008-2014 Johan Hake
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Anders Logg, 2008-2009.
# Modified by Martin Sandve Alnæs 2013-2014

from __future__ import print_function

__all__ = ["Expression"]

# FIXME: Make all error messages uniform according to the following template:
#
# if not isinstance(foo, Foo):
#     raise TypeError("Illegal argument for creation of Bar, not a Foo: " + str(foo))

# Python imports
import types
from six import add_metaclass # Requires newer six version than some buildbots have
from six import string_types
from six.moves import xrange as range
from functools import reduce
import weakref

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
