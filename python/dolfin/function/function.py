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

import ufl
import dolfin.cpp as cpp


class Function(ufl.Coefficient):

    def __init__(self, *args, **kwargs):
        """Initialize Function."""

        if isinstance(args[0], cpp.function.FunctionSpace):
            V = args[0]

            # If initialising from a FunctionSpace
            if len(args) == 1:
                # If passing only the FunctionSpace
                self._cpp_function = cpp.function.Function(V)

                # Initialize the ufl.FunctionSpace
                ufl.Coefficient.__init__(self, V.ufl_function_space(), count=self._cpp_function.id())

        else:
            raise TypeError("expected a FunctionSpace or a Function as argument 1")

    def interpolate(self, u):
        self._cpp_function.interpolate(u)

    def vector(self):
        return self._cpp_function.vector()
