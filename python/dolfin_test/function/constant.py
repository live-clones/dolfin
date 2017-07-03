# -*- coding: utf-8 -*-
"""Create a constant-valued function with given value."""

# Copyright (C) 2008-2015 Anders Logg
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
# Modified by Johan Hake, 2008.
# Modified by Martin Sandve Alnæpps, 2008.

from __future__ import print_function

#__all__ = ["Constant"]

# Import UFL and SWIG-generated extension module (DOLFIN C++)
import ufl
import dolfin_test.cpp as cpp
import numpy

class Constant(ufl.Coefficient, cpp.function.Constant):

    def __init__(self, value, cell=None, name=None):
        """
        Create constant-valued function with given value.

        *Arguments*
            value
                The value may be either a single scalar value, or a
                tuple/list of values for vector-valued functions, or
                nested lists or a numpy array for tensor-valued
                functions.
            cell
                Optional argument. A :py:class:`Cell
                <ufl.Cell>` which defines the geometrical
                dimensions the Constant is defined for.
            name
                Optional argument. A str which overrules the default
                name of the Constant.

        The data type Constant represents a constant value that is
        unknown at compile-time. Its values can thus be changed
        without requiring re-generation and re-compilation of C++
        code.

        *Examples of usage*

            .. code-block:: python

                p = Constant(pi/4)              # scalar
                C = Constant((0.0, -1.0, 0.0))  # constant vector

        """

        # TODO: Either take mesh instead of cell, or drop cell and let
        # grad(c) be undefined.
        if cell is not None:
            cell = ufl.as_cell(cell)
        ufl_domain = None

        array = numpy.array(value)
        rank = len(array.shape)
        floats = list(map(float, array.flat))

        # Create UFL element and initialize constant
        if rank == 0:
            ufl_element = ufl.FiniteElement("Real", cell, 0)
            cpp.function.Constant.__init__(self, floats[0])
        elif rank == 1:
            ufl_element = ufl.VectorElement("Real", cell, 0, dim=len(floats))
            cpp.function.Constant.__init__(self, floats)
        else:
            ufl_element = ufl.TensorElement("Real", cell, 0, shape=array.shape)
            cpp.function.Constant.__init__(self, list(array.shape), floats)

        # Initialize base classes
        ufl_function_space = ufl.FunctionSpace(ufl_domain, ufl_element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())

        # Set name as given or automatic
        name = name or "f_%d" % self.count()
        self.rename(name, "a Constant")

    def cell(self):
        return self.ufl_element().cell()

    def __float__(self):
        # Overriding UFL operator in this particular case.
        if self.ufl_shape:
            raise TypeError("Cannot convert nonscalar constant to float.")
        return cpp.Constant.__float__(self)

    def __str__(self):
        "x.__str__() <==> print(x)"
        return self.name()