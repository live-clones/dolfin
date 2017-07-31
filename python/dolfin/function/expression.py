# -*- coding: utf-8 -*-
from __future__ import print_function

__all__ = ["CompiledExpression", "UserExpression"]

# Python imports
import abc
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
    """A DOLFIN C++ Expression to which user eval functions are attached.

    """

    def __init__(self, user_expression):
        self.user_expression = user_expression

        # Wrap eval functions
        def wrapped_eval(self, values, x):
            self.user_expression.eval(values, x)
        def wrapped_eval_cell(self, values, x, cell):
            self.user_expression.eval_cell(values, x, cell)

        # Attach user-provied Python eval functions (if they exist in
        # the user expression class) to the C++ class
        if hasattr(user_expression, 'eval'):
            print("*** Attaching eval")
            self.eval = types.MethodType(wrapped_eval, self)
        elif hasattr(user_expression, 'eval_cell'):
            print("*** Attaching eval_cell")
            self.eval_cell = types.MethodType(wrapped_eval_cell, self)

        # Create C++ Expression object
        cpp.function.Expression.__init__(self)



class BaseExpression(ufl.Coefficient):
    def __init__(self, cell=None, element=None):

        # Initialise base class
        ufl_function_space = ufl.FunctionSpace(None, element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())

        #name = name or "f_" + str(ufl.Coefficient.count(self))
        #label = label or "User defined expression"
        #self._cpp_object.rename(name, label)

    #@abc.abstractproperty
    #def value(self):
    #    return 'Should never get here'

    def __call__(self, x):
        return self._cpp_object(x)

    def id(self):
        return self._cpp_object.id()

    def value_rank(self):
        return self._cpp_object.value_rank()

    def value_dimension(self, i):
        return self._cpp_object.value_dimension(i)

    def name(self):
        return self._cpp_object.name()

    def label(self):
        return self._cpp_object.label()

    def __str__(self):
        return self._cpp_object.name()

    def cpp_object(self):
        return self._cpp_object

    def compute_vertex_values(self, mesh):
        return self._cpp_object.compute_vertex_values(mesh)


class UserExpression(BaseExpression):
    """Base class for user-defined Python Expression classes, wherer the
    user overload eval or eval_cell

    """

    def __init__(self, *args, **kwargs):

        self._cpp_object = _InterfaceExpression(self)

        # Extract data
        arguments = ("element", "degree", "cell", "domain", "name", "label", "mpi_comm")
        element = kwargs.get("element", None)
        degree = kwargs.get("degree", None)
        cell = kwargs.get("cell", None)
        domain = kwargs.get("domain", None)
        name = kwargs.get("name", None)
        label = kwargs.get("label", None)
        mpi_comm = kwargs.get("mpi_comm", None)

        # Deduce element type if not provided
        if element is None:
            value_shape = tuple(self.value_dimension(i)
                                for i in range(self.value_rank()))
            element = _select_element(family=None, cell=None, degree=2,
                                      value_shape=value_shape)

        BaseExpression.__init__(self, cell=None, element=element)


class Expression(BaseExpression):
    """JIT Expressions"""

    def __init__(self, cpp_code=None, *args, **kwargs):

        # Extract data
        arguments = ("element", "degree", "cell", "domain", "name", "label", "mpi_comm")
        element = kwargs.get("element", None)
        degree = kwargs.get("degree", None)
        cell = kwargs.get("cell", None)
        domain = kwargs.get("domain", None)
        name = kwargs.get("name", None)
        label = kwargs.get("label", None)
        mpi_comm = kwargs.get("mpi_comm", None)

        # Remove arguments that are used in Expression creation
        properties = dict(item for item in kwargs.items() if item[0]  not in arguments)
        for k in properties:
            if not isinstance(k, string_types):
                raise KeyError("Invalid key")
            if not isinstance(properties[k], float):
                raise ValueError("Invalid value")

        self._cpp_object = jit.compile_expression(cpp_code, properties)


        # Deduce element type if not provided
        if element is None:
            value_shape = tuple(self.value_dimension(i)
                                for i in range(self.value_rank()))
            element = _select_element(family=None, cell=None, degree=2,
                                      value_shape=value_shape)


        BaseExpression.__init__(self, cell=None, element=element)

    # This is added dynamically in the intialiser to allow checking of
    # eval in user classes.
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


# Temporary alias for CompiledExpression name
class CompiledExpression(Expression):
    def __init__(self, *args, **kwargs):
        super(CompiledExpression, self).__init__(*args, **kwargs)
