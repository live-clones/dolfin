# -*- coding: utf-8 -*-
"""This module handles the Function class in Python.
"""
# Copyright (C) 2009-2014 Johan Hake
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
# Modified by Martin Sandve Alnæs 2013-2014
# Modified by Anders Logg 2015

#__all__ = ["Function", "TestFunction", "TrialFunction", "Argument",
#           "TestFunctions", "TrialFunctions"]

from six import string_types
import types
import numpy as np

import ufl
import dolfin.cpp as cpp
import dolfin.la as la


class Function(ufl.Coefficient):

    def __init__(self, *args, **kwargs):
        """Initialize Function."""

        if isinstance(args[0], cpp.function.FunctionSpace):
            V = args[0]

            # If initialising from a FunctionSpace
            if len(args) == 1:
                # If passing only the FunctionSpace
                self._cpp_object = cpp.function.Function(V)

                # Initialize the ufl.FunctionSpace
                ufl.Coefficient.__init__(self, V.ufl_function_space(), count=self._cpp_object.id())

        else:
            raise TypeError("expected a FunctionSpace or a Function as argument 1")

    def __call__(self, x):
        return self._cpp_object.__call__(x)

    def interpolate(self, u):
        if isinstance(u, ufl.Coefficient):
            self._cpp_object.interpolate(u._cpp_object)
        else:
            self._cpp_object.interpolate(u)

    def function_space(self):
        return self._cpp_object.function_space()

    def vector(self):
        return self._cpp_object.vector()

    def __float__(self):
        # FIXME: this could be made simple on the C++ (in particular,
        # with dolfin::Scalar)
        if self.ufl_shape != ():
            raise RuntimeError("Cannot convert nonscalar function to float.")
        elm = self.ufl_element()
        if elm.family() != "Real":
            raise RuntimeError("Cannot convert spatially varying function to float.")
        # FIXME: This could be much simpler be exploiting that the
        # vector is ghosted
        # Gather value directly from vector in a parallel safe way
        vec = self.vector()
        indices = np.zeros(1, dtype=la.la_index_dtype())
        values = vec.gather(indices)
        return float(values[0])

    def cpp_object(self):
        return self._cpp_object
