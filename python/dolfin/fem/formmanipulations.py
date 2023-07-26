# -*- coding: utf-8 -*-
# Copyright (C) 2010-2012 Marie E. Rognes
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

import ufl_legacy as ufl
from ufl_legacy.algorithms import elementtransformations, formsplitter

from ufl_legacy.geometry import SpatialCoordinate
from dolfin.function.functionspace import FunctionSpace
from dolfin.function.function import Function
from dolfin.function.argument import Argument


__all__ = ["derivative", "adjoint", "increase_order", "tear", "extract_blocks"]


def adjoint(form, reordered_arguments=None):

    # Call UFL directly if new arguments are provided directly
    if reordered_arguments is not None:
        return ufl.adjoint(form, reordered_arguments=reordered_arguments)

    # Extract form arguments
    arguments = form.arguments()
    if any(arg.part() is not None for arg in arguments):
        raise RuntimeError("Compute adjoint of form, parts not supported")

    if not (len(arguments) == 2):
        raise RuntimeError("Compute adjoint of form, form is not bilinear")

    # Define new Argument(s) in the same spaces (NB: Order does not
    # matter anymore here because number is absolute)
    v_1 = Argument(arguments[1].function_space(), arguments[0].number(),
                   arguments[0].part())
    v_0 = Argument(arguments[0].function_space(), arguments[1].number(),
                   arguments[1].part())

    # Call ufl.adjoint with swapped arguments as new arguments
    return ufl.adjoint(form, reordered_arguments=(v_1, v_0))


adjoint.__doc__ = ufl.adjoint.__doc__


def derivative(form, u, du=None, coefficient_derivatives=None):
    if du is None:
        # Get existing arguments from form and position the new one
        # with the next argument number
        form_arguments = form.arguments()
        number = max([-1] + [arg.number() for arg in form_arguments]) + 1

        # NOTE : Mixed-domains problems need to have arg.part() != None
        # if any(arg.part() is not None for arg in form_arguments):
        #     raise RuntimeError("Compute derivative of form, cannot automatically create new Argument using parts, please supply one")
        part = None

        if isinstance(u, Function):
            # u.part() is None except with mixed-domains
            part = u.part()
            V = u.function_space()
            du = Argument(V, number, part)
        elif isinstance(u, SpatialCoordinate):
            mesh = u.ufl_domain().ufl_cargo()
            element = u.ufl_domain().ufl_coordinate_element()
            V = FunctionSpace(mesh, element)
            du = Argument(V, number, part)
        elif isinstance(u, (list, tuple)) and all(isinstance(w, Function) for w in u):
            raise RuntimeError("Taking derivative of form w.r.t. a tuple of Coefficients. Take derivative w.r.t. a single Coefficient on a mixed space instead.")
        else:
            raise RuntimeError("Computing derivative of form w.r.t. '{}'. Supply Function as a Coefficient".format(u))
    return ufl.derivative(form, u, du, coefficient_derivatives)


derivative.__doc__ = ufl.derivative.__doc__
derivative.__doc__ += """

    A tuple of Coefficients in place of a single Coefficient is not
    supported in DOLFIN. Supply rather a Function on a mixed space in
    place of a Coefficient.
    """


def increase_order(V):
    """For a given function space, return the same space, but with a
    higher polynomial degree

    """
    mesh = V.mesh()
    element = elementtransformations.increase_order(V.ufl_element())
    constrained_domain = V.dofmap().constrained_domain
    return FunctionSpace(mesh, element, constrained_domain=constrained_domain)


def change_regularity(V, family):
    """For a given function space, return the corresponding space with
    the finite elements specified by 'family'. Possible families are
    the families supported by the form compiler

    """
    mesh = V.mesh()
    element = elementtransformations.change_regularity(V.ufl_element(), family)
    constrained_domain = V.dofmap().constrained_domain
    return FunctionSpace(mesh, element, constrained_domain=constrained_domain)


def tear(V):
    """
    For a given function space, return the corresponding discontinuous
    space
    """
    return change_regularity(V, "DG")


def extract_blocks(form, i=None, j=None):
    return formsplitter.extract_blocks(form, i, j)
