# -*- coding: utf-8 -*-
from __future__ import print_function

__all__ = ["CompiledExpression", "UserExpression"]

# Python imports
import hashlib
from functools import reduce
from six import add_metaclass, string_types
from six.moves import xrange as range
import types
import weakref

import dijitso

# Import UFL and SWIG-generated extension module (DOLFIN C++)
import ufl
from ufl import product
from ufl.utils.indexflattening import flatten_multiindex, shape_to_strides
import dolfin.cpp as cpp
import numpy

import dolfin.function.jit as jit

#from dolfin import warning, error

def _select_element(family, cell, degree, value_shape):
    """Select finite element type for cases where user has not provided a
    complete ufl.FiniteElement

    """
    if family is None:
        if degree == 0:
            family = "Discontinuous Lagrange"
        else:
            family = "Lagrange"

    if len(value_shape) == 0:
        element = ufl.FiniteElement(family, cell, degree)
    elif len(value_shape) == 1:
        element = ufl.VectorElement(family, cell, degree, dim=value_shape[0])
    else:
        element = ufl.TensorElement(family, cell, degree, shape=value_shape)

    return element


class _InterfaceExpression(cpp.function.Expression):
    """A DOLFIN C++ Expression . . . . ."""

    def __init__(self, user_expression, *args, **kwargs):
        self.user_expression = user_expression

        cpp.function.Expression.__init__(self, *args, **kwargs)

        def eval(self, values, x):
            self.user_expression.eval(values, x)
        def eval_cell(self, values, x, cell):
            self.user_expression.eval(values, x, cell)

        # Attach eval functions if they exists in the user expression
        # class
        if hasattr(user_expression, 'eval'):
            self.eval = types.MethodType(eval, self)
        if hasattr(user_expression, 'eval_cell'):
            self.eval_cell = types.MethodType(eval_cell, self)


class UserExpression(ufl.Coefficient):
    def __init__(self, function_space=None, element=None, degree=None):
        """Create an Expression."""

        #if element is None and degree is None:
        #    raise RuntimeError('UserExpression must specific a FiniteElement or a dgeree')

        #if element is None:
        #    # Create UFL element
        #    element = ufl.FiniteElement(family, mesh.ufl_cell(), degree,
        #                                form_degree=None)

        #ufl.Coefficient.__init__(self, function_space)
        #cpp.function.Expression.__init__(self, 1)
        #cpp.function.Expression.__init__(self, self.ufl_shape)

        self._cpp_object = _InterfaceExpression(self)
        value_shape = tuple(self.value_dimension(i)
                            for i in range(self.value_rank()))
        if element is None:
            element = _select_element(family=None, cell=None, degree=2,
                                      value_shape=value_shape)

        # Initialize UFL base class
        ufl_function_space = ufl.FunctionSpace(None, element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())

    def value_rank(self):
        return self._cpp_object.value_rank()

    def value_dimension(self, i):
        return self._cpp_object.value_dimension(i)

    def id(self):
        return self._cpp_object.id()

    def cpp_object(self):
        """Return the underling cpp.Expression object"""
        return self._cpp_object


class CompiledExpression(ufl.Coefficient):
    def __init__(self, statements, **kwargs):

        # Extract data
        degree = kwargs.pop("degree", None)
        element = kwargs.pop("element", None)

        # Determine Expression type (JIT or user overlaoded)

        # Deduce underlying element is not explicitly provided

        properties = kwargs
        for k in properties:
            if not isinstance(k, string_types):
                raise KeyError("Invalid key")
            if not isinstance(properties[k], float):
                raise ValueError("Invalid value")

        self._cpp_object = jit.compile_expression(statements, properties)

        if element is None:
            value_shape = tuple(self.value_dimension(i)
                                for i in range(self.value_rank()))
            element = _select_element(family=None, cell=None, degree=2,
                                      value_shape=value_shape)

        ufl_function_space = ufl.FunctionSpace(None, element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())

    def __getattr__(self, name):
        "Pass attributes through to (JIT compiled) Expression object"
        if hasattr(self._cpp_object, name):
            return self._cpp_object.get_property(name)
        else:
            raise(AttributeError)

    def __setattr__(self, name, value):
        if name.startswith("_"):
            super().__setattr__(name, value)
        else:
            self._cpp_object.set_property(name, value)

    def __call__(self, x):
        return self._cpp_object(x)

    def id(self):
        return self._cpp_object.id()

    def value_rank(self):
        return self._cpp_object.value_rank()

    def value_dimension(self, i):
        return self._cpp_object.value_dimension(i)

    def cpp_object(self):
        return self._cpp_object

    def compute_vertex_values(self, mesh):
        return self._cpp_object.compute_vertex_values(mesh)
