# -*- coding: utf-8 -*-

# Copyright (C) 2011-2017 Anders Logg and Garth N. Wells
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

from ufl.form import sub_forms_by_domain

import dolfin.cpp as cpp
import dolfin.fem.solving
from dolfin.fem.form import Form


class LinearVariationalProblem(cpp.fem.LinearVariationalProblem):

    def __init__(self, a, L, u, bcs=None, form_compiler_parameters=None):
        """Create linear variational problem a(u, v) = L(v).

        An optional argument bcs may be passed to specify boundary
        conditions.

        Another optional argument form_compiler_parameters may be
        specified to pass parameters to the form compiler.

        """

        if bcs is None:
            bcs = []
        elif not isinstance(bcs, (list, tuple)):
            bcs = [bcs]

        # Store input UFL forms and solution Function
        self.a_ufl = a
        self.L_ufl = L
        self.u_ufl = u

        # Store form compiler parameters
        form_compiler_parameters = form_compiler_parameters or {}
        self.form_compiler_parameters = form_compiler_parameters

        # Wrap forms (and check if linear form L is empty)
        if L.empty():
            L = cpp.fem.Form(1, 0)
        else:
            L = Form(L, form_compiler_parameters=form_compiler_parameters)
        a = Form(a, form_compiler_parameters=form_compiler_parameters)

        # Initialize C++ base class
        cpp.fem.LinearVariationalProblem.__init__(self, a, L, u._cpp_object, bcs)


class MixedLinearVariationalProblem(cpp.fem.MixedLinearVariationalProblem):

    def __init__(self, a, L, u, bcs=None, form_compiler_parameters=None):
        """Create mixed linear variational problem a(u, v) = L(v).

        An optional argument bcs may be passed to specify boundary
        conditions.

        Another optional argument form_compiler_parameters may be
        specified to pass parameters to the form compiler.

        """

        # Extract and check arguments (u is a list of Function)
        u_comps = [u[i]._cpp_object for i in range(len(u))]
        bcs = dolfin.fem.solving._extract_bcs(bcs)

        # Store form compiler parameters
        form_compiler_parameters = form_compiler_parameters or {}
        self.form_compiler_parameters = form_compiler_parameters

        # Update rhs if we don't have a consistent number of blocks
        if len(L) != len(u):
            L_tmp = [None for i in range(len(u))]
            for Li in L:
                L_tmp[Li.arguments()[0].part()] = Li
            L = L_tmp

        # Check number of blocks in lhs, rhs are consistent
        assert(len(a) == len(u) * len(u))
        assert(len(L) == len(u))

        # Create list of forms/blocks
        a_list = list()
        L_list = list()
        for Li in L:
            if Li is None:
                L_list.append([cpp.fem.Form(1, 0)])  # single-elt list
            elif Li.empty():
                L_list.append([cpp.fem.Form(1, 0)])  # single-elt list
            else:
                Ls = []  # List of Li subforms
                for Lsub in sub_forms_by_domain(Li):
                    if Lsub is None:
                        Ls.append(cpp.fem.Form(1, 0))
                    elif Lsub.empty():
                        Ls.append(cpp.fem.Form(1, 0))
                    else:
                        Ls.append(Form(Lsub, form_compiler_parameters=form_compiler_parameters))
                L_list.append(Ls)

        for ai in a:
            if ai is None:
                a_list.append([cpp.fem.Form(2, 0)])
            else:
                As = []
                for Asub in sub_forms_by_domain(ai):
                    As.append(Form(Asub, form_compiler_parameters=form_compiler_parameters))
                a_list.append(As)

        # Initialize C++ base class
        cpp.fem.MixedLinearVariationalProblem.__init__(self, a_list, L_list, u_comps, bcs)


class NonlinearVariationalProblem(cpp.fem.NonlinearVariationalProblem):

    def __init__(self, F, u, bcs=None, J=None, form_compiler_parameters=None):
        """Create nonlinear variational problem F(u; v) = 0.

        Optional arguments bcs and J may be passed to specify boundary
        conditions and the Jacobian J = dF/du.

        Another optional argument form_compiler_parameters may be
        specified to pass parameters to the form compiler.

        """

        if bcs is None:
            bcs = []
        elif not isinstance(bcs, (list, tuple)):
            bcs = [bcs]

        # Store input UFL forms and solution Function
        self.F_ufl = F
        self.J_ufl = J
        self.u_ufl = u

        # Store form compiler parameters
        form_compiler_parameters = form_compiler_parameters or {}
        self.form_compiler_parameters = form_compiler_parameters

        # Wrap forms
        F = Form(F, form_compiler_parameters=form_compiler_parameters)
        if J is not None:
            J = Form(J, form_compiler_parameters=form_compiler_parameters)

        # Initialize C++ base class
        cpp.fem.NonlinearVariationalProblem.__init__(self, F, u._cpp_object, bcs, J)
