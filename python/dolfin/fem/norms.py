# -*- coding: utf-8 -*-
"""This module provides a simple way to compute various norms of
:py:class:`Vectors <dolfin.cpp.Vector>` and :py:class:`Functions
<dolfin.functions.function.Function>`, including the standard
:math:`L^2`-norm and other norms.

"""

# Copyright (C) 2008-2014 Anders Logg
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

from math import sqrt
import ufl_legacy as ufl
import functools
from ufl_legacy import (grad, div, curl)
import dolfin.cpp as cpp
from dolfin.fem.assembling import assemble, assemble_multimesh
from dolfin.fem.interpolation import interpolate
from dolfin.function.functionspace import (FunctionSpace,
                                           VectorFunctionSpace,
                                           TensorFunctionSpace)
from dolfin.function.function import Function
from dolfin.function.multimeshfunction import MultiMeshFunction
from dolfin.function.multimeshfunctionspace import (MultiMeshFunctionSpace,
                                                    MultiMeshVectorFunctionSpace, MultiMeshTensorFunctionSpace)

__all__ = ["norm", "errornorm"]


def norm(v, norm_type="L2", mesh=None):
    r"""
    Return the norm of a given vector or function.

    *Arguments*
        v
            a :py:class:`Vector <dolfin.cpp.Vector>` or
            a :py:class:`Function <dolfin.functions.function.Function>`.
        norm_type
            see below for alternatives.
        mesh
            optional :py:class:`Mesh <dolfin.cpp.Mesh>` on
            which to compute the norm.

    If the norm type is not specified, the standard :math:`L^2` -norm
    is computed. Possible norm types include:

    *Vectors*

    ================   =================  ================
    Norm               Usage
    ================   =================  ================
    :math:`l^2`        norm(x, 'l2')      Default
    :math:`l^1`        norm(x, 'l1')
    :math:`l^\infty`   norm(x, 'linf')
    ================   =================  ================

    *Functions*

    ================  =================  =================================
    Norm              Usage              Includes the :math:`L^2` -term
    ================  =================  =================================
    :math:`L^2`       norm(v, 'L2')      Yes
    :math:`H^1`       norm(v, 'H1')      Yes
    :math:`H^1_0`     norm(v, 'H10')     No
    :math:`H` (div)   norm(v, 'Hdiv')    Yes
    :math:`H` (div)   norm(v, 'Hdiv0')   No
    :math:`H` (curl)  norm(v, 'Hcurl')   Yes
    :math:`H` (curl)  norm(v, 'Hcurl0')  No
    ================  =================  =================================

    *Examples of usage*

    .. code-block:: python

        v = Function(V)
        x = v.vector()

        print norm(x, 'linf')   # print the infinity norm of vector x

        n = norm(v)             # compute L^2 norm of v
        print norm(v, 'Hdiv')   # print H(div) norm of v
        n = norm(v, 'H1', mesh) # compute H^1 norm of v on given mesh

    """

    # if not isinstance(v, (GenericVector, GenericFunction)):
    #     cpp.dolfin_error("norms.py",
    #                      "compute norm",
    #                      "expected a GenericVector or GenericFunction")

    # Check arguments
    # if not isinstance(norm_type, string_types):
    #     cpp.dolfin_error("norms.py",
    #                      "compute norm",
    #                      "Norm type must be a string, not " +
    #                      str(type(norm_type)))
    # if mesh is not None and not isinstance(mesh, cpp.Mesh):
    #     cpp.dolfin_error("norms.py",
    #                      "compute norm",
    #                      "Expecting a Mesh, not " + str(type(mesh)))

    # Get mesh from function
    if isinstance(v, cpp.function.Function) and mesh is None:
        mesh = v.function_space().mesh()
    elif isinstance(v, MultiMeshFunction) and mesh is None:
        mesh = v.function_space().multimesh()

    # Define integration measure and domain
    if isinstance(v, MultiMeshFunction):
        dc = ufl.dx(mesh) + ufl.dC(mesh)
        assemble_func = functools.partial(assemble_multimesh, form_compiler_parameters={"representation": "quadrature"})
    else:
        dc = ufl.dx(mesh)
        assemble_func = assemble
    # Select norm type
    if isinstance(v, cpp.la.GenericVector):
        return v.norm(norm_type.lower())

    elif isinstance(v, ufl.Coefficient):
        if norm_type.lower() == "l2":
            M = v**2 * dc
        elif norm_type.lower() == "h1":
            M = (v**2 + grad(v)**2) * dc
        elif norm_type.lower() == "h10":
            M = grad(v)**2 * dc
        elif norm_type.lower() == "hdiv":
            M = (v**2 + div(v)**2) * dc
        elif norm_type.lower() == "hdiv0":
            M = div(v)**2 * dc
        elif norm_type.lower() == "hcurl":
            M = (v**2 + curl(v)**2) * dc
        elif norm_type.lower() == "hcurl0":
            M = curl(v)**2 * dc
        else:
            raise ValueError("Unknown norm type {}".format(str(norm_type)))
    else:
        raise TypeError("Do not know how to compute norm of {}".format(str(v)))

    # Assemble value and return
    return sqrt(assemble_func(M))


def errornorm(u, uh, norm_type="l2", degree_rise=3, mesh=None):
    """
    Compute and return the error :math:`e = u - u_h` in the given norm.

    *Arguments*
        u, uh
            :py:class:`Functions <dolfin.functions.function.Function>`
        norm_type
            Type of norm. The :math:`L^2` -norm is default.
            For other norms, see :py:func:`norm <dolfin.fem.norms.norm>`.
        degree_rise
            The number of degrees above that of u_h used in the
            interpolation; i.e. the degree of piecewise polynomials used
            to approximate :math:`u` and :math:`u_h` will be the degree
            of :math:`u_h` + degree_raise.
        mesh
            Optional :py:class:`Mesh <dolfin.cpp.Mesh>` on
            which to compute the error norm.

    In simple cases, one may just define

    .. code-block:: python

        e = u - uh

    and evalute for example the square of the error in the :math:`L^2` -norm by

    .. code-block:: python

        assemble(e**2*dx(mesh))

    However, this is not stable w.r.t. round-off errors considering that
    the form compiler may expand(#) the expression above to::

        e**2*dx = u**2*dx - 2*u*uh*dx + uh**2*dx

    and this might get further expanded into thousands of terms for
    higher order elements. Thus, the error will be evaluated by adding
    a large number of terms which should sum up to something close to
    zero (if the error is small).

    This module computes the error by first interpolating both
    :math:`u` and :math:`u_h` to a common space (of high accuracy),
    then subtracting the two fields (which is easy since they are
    expressed in the same basis) and then evaluating the integral.

    (#) If using the tensor representation optimizations.
    The quadrature represenation does not suffer from this problem.
    """

    # Check argument
    # if not isinstance(u, cpp.function.GenericFunction):
    #     cpp.dolfin_error("norms.py",
    #                      "compute error norm",
    #                      "Expecting a Function or Expression for u")
    # if not isinstance(uh, cpp.function.Function):
    #     cpp.dolfin_error("norms.py",
    #                      "compute error norm",
    #                      "Expecting a Function for uh")

    # Get mesh
    if isinstance(u, cpp.function.Function) and mesh is None:
        mesh = u.function_space().mesh()
    if isinstance(uh, cpp.function.Function) and mesh is None:
        mesh = uh.function_space().mesh()
    if isinstance(uh, MultiMeshFunction) and mesh is None:
        mesh = uh.function_space().multimesh()
    if hasattr(uh, "_cpp_object") and mesh is None:
        mesh = uh._cpp_object.function_space().mesh()
    if hasattr(u, "_cpp_object") and mesh is None:
        mesh = u._cpp_object.function_space().mesh()
    if mesh is None:
        raise RuntimeError("Cannot compute error norm. Missing mesh.")

    # Get rank
    if not u.ufl_shape == uh.ufl_shape:
        raise RuntimeError("Cannot compute error norm. Value shapes do not match.")

    shape = u.ufl_shape
    rank = len(shape)

    # Check that uh is associated with a finite element
    if uh.ufl_element().degree() is None:
        raise RuntimeError("Cannot compute error norm. Function uh must have a finite element.")

    # Degree for interpolation space. Raise degree with respect to uh.
    degree = uh.ufl_element().degree() + degree_rise

    # Check degree of 'exact' solution u
    degree_u = u.ufl_element().degree()
    if degree_u is not None and degree_u < degree:
        cpp.warning("Degree of exact solution may be inadequate for accurate result in errornorm.")

    # Create function space
    if isinstance(uh, MultiMeshFunction):
        function_space = MultiMeshFunctionSpace
        vector_space = MultiMeshVectorFunctionSpace
        tensor_space = MultiMeshTensorFunctionSpace
        func = MultiMeshFunction
    else:
        function_space = FunctionSpace
        vector_space = VectorFunctionSpace
        tensor_space = TensorFunctionSpace
        func = Function

    if rank == 0:
        V = function_space(mesh, "Discontinuous Lagrange", degree)
    elif rank == 1:
        V = vector_space(mesh, "Discontinuous Lagrange", degree, dim=shape[0])
    elif rank > 1:
        V = tensor_space(mesh, "Discontinuous Lagrange", degree, shape=shape)

    # Interpolate functions into finite element space
    pi_u = interpolate(u, V)
    pi_uh = interpolate(uh, V)

    # Compute the difference
    e = func(V)

    e.assign(pi_u)
    e.vector().axpy(-1.0, pi_uh.vector())

    # Compute norm
    return norm(e, norm_type=norm_type, mesh=mesh)
