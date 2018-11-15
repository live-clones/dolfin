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


class MixedNonlinearVariationalProblem(cpp.fem.MixedNonlinearVariationalProblem):

    def __init__(self, F, u, bcs=None, J=None, form_compiler_parameters=None):
        """Create nonlinear variational problem F(u; v) = 0.

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

        # Store input UFL forms and solution Function
        self.F_ufl = F
        self.J_ufl = J
        self.u_ufl = u

        # Update rhs if we don't have a consistent number of blocks
        if len(F) != len(u):
            F_tmp = [None for i in range(len(u))]
            for Fi in F:
                F_tmp[Fi.arguments()[0].part()] = Fi
            F = F_tmp

        # Check number of blocks in the residual and solution are coherent
        assert(len(J) == len(u) * len(u))
        assert(len(F) == len(u))

        # Create list of forms/blocks
        F_list = list()
        print("[problem.py] size F = ", len(F))
        for Fi in F:
            if Fi is None:
                F_list.append([cpp.fem.Form(1, 0)])
            elif Fi.empty():
                F_list.append([cpp.fem.Form(1, 0)])  # single-elt list
            else:
                Fs = []
                for Fsub in sub_forms_by_domain(Fi):
                    if Fsub is None:
                        Fs.append(cpp.fem.Form(1, 0))
                    elif Fsub.empty():
                        Fs.append(cpp.fem.Form(1, 0))
                    else:
                        Fs.append(Form(Fsub, form_compiler_parameters=form_compiler_parameters))
                F_list.append(Fs)
        print("[problem] create list of residual forms OK")

        J_list = None
        if J is not None:
            J_list = list()
            print("[problem.py] size J = ", len(J))
            for Ji in J:
                if Ji is None:
                    J_list.append([cpp.fem.Form(2, 0)])
                elif Ji.empty():
                    J_list.append([cpp.fem.Form(2, 0)])
                else:
                    Js = []
                    for Jsub in sub_forms_by_domain(Ji):
                        Js.append(Form(Jsub, form_compiler_parameters=form_compiler_parameters))
                    J_list.append(Js)
        print("[problem] create list of jacobian forms OK, J_list size = ", len(J_list))

        # Initialize C++ base class
        cpp.fem.MixedNonlinearVariationalProblem.__init__(self, F_list, u_comps, bcs, J_list)
