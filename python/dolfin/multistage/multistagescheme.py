# -*- coding: utf-8 -*-
"""This module defines different MultiStageScheme classes which can be
passed to a RKSolver or PointIntegralSolver

"""

# Copyright (C) 2013-2015 Johan Hake
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
# Modified by Patrick Farrell, 2013
# Modified by Martin Sandve Alnæs, 2015

import numpy as np
import functools
import ufl_legacy as ufl

import dolfin.cpp as cpp
from dolfin.function.constant import Constant
from dolfin.function.expression import Expression
from dolfin.function.function import Function
from dolfin.fem.formmanipulations import derivative, adjoint
from ufl_legacy import action as ufl_action
from dolfin.fem.form import Form
import ufl_legacy.algorithms
from ufl_legacy.algorithms import expand_derivatives

# FIXME: Add support for algebraic parts (at least for implicit)
# FIXME: Add support for implicit/explicit split ala IMEX schemes


def safe_adjoint(x):
    return adjoint(x, reordered_arguments=x.arguments())


def safe_action(x, y):
    x = expand_derivatives(x)
    if x.integrals() == ():
        return x  # form is empty, return anyway
    else:
        return ufl_action(x, y)


def _check_abc(a, b, c):
    if not (isinstance(a, np.ndarray) and (len(a) == 1 or
            (len(a.shape) == 2 and a.shape[0] == a.shape[1]))):
        raise TypeError("Expected an m x m numpy array as the first argument")
    if not (isinstance(b, np.ndarray) and len(b.shape) in [1, 2]):
        raise TypeError("Expected a 1 or 2 dimensional numpy array as the second argument")
    if not (isinstance(c, np.ndarray) and len(c.shape) == 1):
        raise TypeError("Expected a 1 dimensional numpy array as the third argument")

    # Make sure a is a "matrix"
    if len(a) == 1:
        a.shape = (1, 1)

    # Get size of system
    size = a.shape[0]

    # If b is a matrix we expect it to have two rows
    if len(b.shape) == 2:
        if not (b.shape[0] == 2 and b.shape[1] == size):
            raise ValueError("Expected a 2 row matrix with the same number "
                             "of collumns as the first dimension of the a matrix.")
    elif len(b) != size:
        raise ValueError("Expected the length of the b vector to have the "
                         "same size as the first dimension of the a matrix.")

    if len(c) != size:
        raise ValueError("Expected the length of the c vector to have the "
                         "same size as the first dimension of the a matrix.")

    # Check if the method is singly diagonally implicit
    sigma = -1
    for i in range(size):
        # If implicit
        if a[i, i] != 0:
            if sigma == -1:
                sigma = a[i, i]
            elif sigma != a[i, i]:
                raise ValueError("Expected only singly diagonally implicit "
                                 "schemes. (Same value on the diagonal of 'a'.)")

    # Check if tableau is fully implicit
    for i in range(size):
        for j in range(i):
            if a[j, i] != 0:
                raise ValueError("Does not support fully implicit Butcher tableau.")

    return a


def _check_form(rhs_form):
    if not isinstance(rhs_form, ufl.Form):
        raise TypeError("Expected a ufl.Form as the 5th argument.")

    # Check if form contains a cell or point integral
    if rhs_form.integrals_by_type("cell"):
        DX = ufl.dx
    elif rhs_form.integrals_by_type("vertex"):
        DX = ufl.dP
    else:
        raise ValueError("Expected either a cell or vertex integral in the form.")

    if len(rhs_form.integrals()) != 1:
        raise ValueError("Expected only one integral in form.")

    arguments = rhs_form.arguments()
    if len(arguments) != 1:
        raise ValueError("Expected the form to have rank 1")

    return DX


def _time_dependent_expressions(rhs_form, time):
    """Return a list of expressions which uses the present time as a
    parameter

    """
    # FIXME: Add extraction of time dependant expressions from bcs too
    time_dependent_expressions = dict()

    for coefficient in rhs_form.coefficients():
        if hasattr(coefficient, "_user_parameters"):
            for c_name, c in list(coefficient._user_parameters.items()):
                if isinstance(c, ufl.Coefficient) and time.id() == c.id():
                    if coefficient not in time_dependent_expressions:
                        time_dependent_expressions[coefficient] = [c_name]
                    else:
                        time_dependent_expressions[coefficient].append(c_name)

    return time_dependent_expressions


def _replace_dict_time_dependent_expression(time_dep_expressions, time,
                                            dt, c):
    assert (isinstance(c, float))
    replace_dict = {}
    if c == 0.0 or not time_dep_expressions:
        return replace_dict
    new_time = Expression("time + c*dt", time=time, c=c, dt=dt, degree=0)
    for expr, c_names in list(time_dep_expressions.items()):
        assert (isinstance(expr, Expression))
        kwargs = dict(name=expr.name(), label=expr.label(),
                      element=expr.ufl_element(), **expr._user_parameters)
        for c_name in c_names:
            kwargs[c_name] = new_time
        replace_dict[expr] = Expression(expr._cppcode, **kwargs)

    return replace_dict


def _butcher_scheme_generator(a, b, c, time, solution, rhs_form):
    """Generates a list of forms and solutions for a given Butcher tableau

    *Arguments*
        a (2 dimensional numpy array)
            The a matrix of the Butcher tableau.
        b (1-2 dimensional numpy array)
            The b vector of the Butcher tableau. If b is 2 dimensional the
            scheme includes an error estimator and can be used in adaptive
            solvers.
        c (1 dimensional numpy array)
            The c vector the Butcher tableau.
        time (_Constant_)
            A Constant holding the time at the start of the time step
        solution (_Function_)
            The prognostic variable
        rhs_form (ufl.Form)
            A UFL form representing the rhs for a time differentiated equation

    """

    a = _check_abc(a, b, c)
    size = a.shape[0]

    DX = _check_form(rhs_form)

    # Get test function
    arguments = rhs_form.arguments()
    # coefficients = rhs_form.coefficients()
    v = arguments[0]

    # Create time step
    dt = Constant(0.1)

    # rhs forms
    dolfin_stage_forms = []
    ufl_stage_forms = []

    # Stage solutions
    k = [Function(solution.function_space(), name="k_%d" % i) for i in range(size)]

    jacobian_indices = []

    # Create the stage forms
    y_ = solution
    time_ = time
    time_dep_expressions = _time_dependent_expressions(rhs_form, time)
    zero_ = ufl.zero(*y_.ufl_shape)
    for i, ki in enumerate(k):

        # Check whether the stage is explicit
        explicit = a[i, i] == 0

        # Evaluation arguments for the ith stage
        evalargs = y_ + dt * sum([float(a[i, j]) * k[j]
                                  for j in range(i + 1)], zero_)
        time = time_ + dt * c[i]

        replace_dict = _replace_dict_time_dependent_expression(time_dep_expressions,
                                                               time_, dt, c[i])

        replace_dict[y_] = evalargs
        replace_dict[time_] = time
        stage_form = ufl.replace(rhs_form, replace_dict)

        if explicit:
            stage_forms = [stage_form]
            jacobian_indices.append(-1)
        else:
            # Create a F=0 form and differentiate it
            stage_form -= ufl.inner(ki, v) * DX
            stage_forms = [stage_form, derivative(stage_form, ki)]
            jacobian_indices.append(0)
        ufl_stage_forms.append(stage_forms)

        dolfin_stage_forms.append([Form(form) for form in stage_forms])

    # Only one last stage
    if len(b.shape) == 1:
        last_stage = Form(ufl.inner(y_ + sum([dt * float(bi) * ki for bi, ki in
                                              zip(b, k)], zero_), v) * DX)
    else:
        # FIXME: Add support for adaptivity in RKSolver and
        # MultiStageScheme
        last_stage = [Form(ufl.inner(y_ + sum([dt * float(bi) * ki for bi, ki in
                                               zip(b[0, :], k)], zero_), v) * DX),
                      Form(ufl.inner(y_ + sum([dt * float(bi) * ki for bi, ki in
                                               zip(b[1, :], k)], zero_), v) * DX)]

    # Create the Function holding the solution at end of time step
    # k.append(solution.copy())

    # Generate human form of MultiStageScheme
    human_form = []
    for i in range(size):
        kterm = " + ".join("%sh*k_%s" % ("" if a[i, j] == 1.0 else
                                         "%s*" % a[i, j], j)
                           for j in range(size) if a[i, j] != 0)
        if c[i] in [0.0, 1.0]:
            cih = " + h" if c[i] == 1.0 else ""
        else:
            cih = " + %s*h" % c[i]

        if len(kterm) == 0:
            human_form.append("k_%(i)s = f(t_n%(cih)s, y_n)" % {"i": i, "cih": cih})
        else:
            human_form.append("k_%(i)s = f(t_n%(cih)s, y_n + %(kterm)s)" %
                              {"i": i, "cih": cih, "kterm": kterm})

    parentheses = "(%s)" if np.sum(b > 0) > 1 else "%s"
    human_form.append("y_{n+1} = y_n + h*" + parentheses % (" + ".join(
        "%sk_%s" % ("" if b[i] == 1.0 else "%s*" % b[i], i)
        for i in range(size) if b[i] > 0)))

    human_form = "\n".join(human_form)

    return ufl_stage_forms, dolfin_stage_forms, jacobian_indices, last_stage, \
        k, dt, human_form, None


def _butcher_scheme_generator_tlm(a, b, c, time, solution, rhs_form,
                                  perturbation):
    """Generates a list of forms and solutions for a given Butcher tableau

    *Arguments*
        a (2 dimensional numpy array)
            The a matrix of the Butcher tableau.
        b (1-2 dimensional numpy array)
            The b vector of the Butcher tableau. If b is 2 dimensional the
            scheme includes an error estimator and can be used in adaptive
            solvers.
        c (1 dimensional numpy array)
            The c vector the Butcher tableau.
        time (_Constant_)
            A Constant holding the time at the start of the time step
        solution (_Function_)
            The prognostic variable
        rhs_form (ufl.Form)
            A UFL form representing the rhs for a time differentiated equation
        perturbation (_Function_)
            The perturbation in the initial condition of the solution

    """

    a = _check_abc(a, b, c)
    size = a.shape[0]

    DX = _check_form(rhs_form)

    # Get test function
    arguments = rhs_form.arguments()
    # coefficients = rhs_form.coefficients()
    v = arguments[0]

    # Create time step
    dt = Constant(0.1)

    # rhs forms
    dolfin_stage_forms = []
    ufl_stage_forms = []

    # Stage solutions
    k = [Function(solution.function_space(), name="k_%d" % i) for i in range(size)]
    kdot = [Function(solution.function_space(), name="kdot_%d" % i)
            for i in range(size)]

    # Create the stage forms
    y_ = solution
    time_ = time
    time_dep_expressions = _time_dependent_expressions(rhs_form, time)
    zero_ = ufl.zero(*y_.ufl_shape)
    forward_forms = []
    stage_solutions = []
    jacobian_indices = []

    for i, ki in enumerate(k):

        # Check whether the stage is explicit
        explicit = a[i, i] == 0

        # Evaluation arguments for the ith stage
        evalargs = y_ + dt * sum([float(a[i, j]) * k[j]
                                  for j in range(i + 1)], zero_)
        time = time_ + dt * c[i]

        replace_dict = _replace_dict_time_dependent_expression(time_dep_expressions,
                                                               time_, dt, c[i])

        replace_dict[y_] = evalargs
        replace_dict[time_] = time
        stage_form = ufl.replace(rhs_form, replace_dict)

        forward_forms.append(stage_form)

        # The recomputation of the forward run:

        if explicit:
            stage_forms = [stage_form]
            jacobian_indices.append(-1)
        else:
            # Create a F=0 form and differentiate it
            stage_form_implicit = stage_form - ufl.inner(ki, v) * DX
            stage_forms = [stage_form_implicit, derivative(stage_form_implicit, ki)]
            jacobian_indices.append(0)

        ufl_stage_forms.append(stage_forms)
        dolfin_stage_forms.append([Form(form) for form in stage_forms])
        stage_solutions.append(ki)

        # And now the tangent linearisation:
        stage_form_tlm = safe_action(derivative(stage_form, y_), perturbation) + \
            sum([dt * float(a[i, j]) * safe_action(derivative(
                forward_forms[j], y_), kdot[j]) for j in range(i + 1)])
        if explicit:
            stage_forms_tlm = [stage_form_tlm]
            jacobian_indices.append(-1)
        else:
            # Create a F=0 form and differentiate it
            stage_form_tlm -= ufl.inner(kdot[i], v) * DX
            stage_forms_tlm = [stage_form_tlm, derivative(stage_form_tlm, kdot[i])]
            jacobian_indices.append(1)

        ufl_stage_forms.append(stage_forms_tlm)
        dolfin_stage_forms.append([Form(form) for form in stage_forms_tlm])
        stage_solutions.append(kdot[i])

    # Only one last stage
    if len(b.shape) == 1:
        last_stage = Form(ufl.inner(perturbation + sum(
            [dt * float(bi) * kdoti for bi, kdoti in zip(b, kdot)], zero_), v) * DX)
    else:
        raise Exception("Not sure what to do here")

    human_form = []
    for i in range(size):
        kterm = " + ".join("%sh*k_%s" % ("" if a[i, j] == 1.0 else
                                         "%s*" % a[i, j], j)
                           for j in range(size) if a[i, j] != 0)
        if c[i] in [0.0, 1.0]:
            cih = " + h" if c[i] == 1.0 else ""
        else:
            cih = " + %s*h" % c[i]

        kdotterm = " + ".join("%(a)sh*action(derivative(f(t_n%(cih)s, y_n + "
                              "%(kterm)s), kdot_%(i)s" %
                              {"a": ("" if a[i, j] == 1.0 else "%s*" % a[i, j], j),
                               "i": i,
                               "cih": cih,
                               "kterm": kterm}
                              for j in range(size) if a[i, j] != 0)

        if len(kterm) == 0:
            human_form.append("k_%(i)s = f(t_n%(cih)s, y_n)" % {"i": i, "cih": cih})
            human_form.append("kdot_%(i)s = action(derivative("
                              "f(t_n%(cih)s, y_n), y_n), ydot_n)" %
                              {"i": i, "cih": cih})
        else:
            human_form.append("k_%(i)s = f(t_n%(cih)s, y_n + %(kterm)s)" %
                              {"i": i, "cih": cih, "kterm": kterm})
            human_form.append("kdot_%(i)s = action(derivative(f(t_n%(cih)s, "
                              "y_n + %(kterm)s), y_n) + %(kdotterm)s" %
                              {"i": i, "cih": cih, "kterm": kterm, "kdotterm": kdotterm})

    parentheses = "(%s)" if np.sum(b > 0) > 1 else "%s"
    human_form.append("ydot_{n+1} = ydot_n + h*" + parentheses % (" + ".join(
        "%skdot_%s" % ("" if b[i] == 1.0 else "%s*" % b[i], i)
        for i in range(size) if b[i] > 0)))

    human_form = "\n".join(human_form)

    return ufl_stage_forms, dolfin_stage_forms, jacobian_indices, last_stage, \
        stage_solutions, dt, human_form, perturbation


def _butcher_scheme_generator_adm(a, b, c, time, solution, rhs_form, adj):
    """Generates a list of forms and solutions for a given Butcher tableau

    *Arguments*
        a (2 dimensional numpy array)
            The a matrix of the Butcher tableau.
        b (1-2 dimensional numpy array)
            The b vector of the Butcher tableau. If b is 2 dimensional the
            scheme includes an error estimator and can be used in adaptive
            solvers.
        c (1 dimensional numpy array)
            The c vector the Butcher tableau.
        time (_Constant_)
            A Constant holding the time at the start of the time step
        solution (_Function_)
            The prognostic variable
        rhs_form (ufl.Form)
            A UFL form representing the rhs for a time differentiated equation
        adj (_Function_)
            The derivative of the functional with respect to y_n+1

    """

    a = _check_abc(a, b, c)
    size = a.shape[0]

    DX = _check_form(rhs_form)

    # Get test function
    arguments = rhs_form.arguments()
    # coefficients = rhs_form.coefficients()
    v = arguments[0]

    # Create time step
    dt = Constant(0.1)

    # rhs forms
    dolfin_stage_forms = []
    ufl_stage_forms = []

    # Stage solutions
    k = [Function(solution.function_space(), name="k_%d" % i) for i in range(size)]
    kbar = [Function(solution.function_space(), name="kbar_%d" % i)
            for i in range(size)]

    # Create the stage forms
    y_ = solution
    time_ = time
    time_dep_expressions = _time_dependent_expressions(rhs_form, time)
    zero_ = ufl.zero(*y_.ufl_shape)
    forward_forms = []
    stage_solutions = []
    jacobian_indices = []

    # The recomputation of the forward run:
    for i, ki in enumerate(k):

        # Check whether the stage is explicit
        explicit = a[i, i] == 0

        # Evaluation arguments for the ith stage
        evalargs = y_ + dt * sum([float(a[i, j]) * k[j]
                                  for j in range(i + 1)], zero_)
        time = time_ + dt * c[i]

        replace_dict = _replace_dict_time_dependent_expression(
            time_dep_expressions, time_, dt, c[i])

        replace_dict[y_] = evalargs
        replace_dict[time_] = time
        stage_form = ufl.replace(rhs_form, replace_dict)

        forward_forms.append(stage_form)

        if explicit:
            stage_forms = [stage_form]
            jacobian_indices.append(-1)
        else:
            # Create a F=0 form and differentiate it
            stage_form_implicit = stage_form - ufl.inner(ki, v) * DX
            stage_forms = [stage_form_implicit, derivative(
                stage_form_implicit, ki)]
            jacobian_indices.append(0)

        ufl_stage_forms.append(stage_forms)
        dolfin_stage_forms.append([Form(form) for form in stage_forms])
        stage_solutions.append(ki)

    for i, kbari in reversed(list(enumerate(kbar))):

        # Check whether the stage is explicit
        explicit = a[i, i] == 0

        # And now the adjoint linearisation:
        stage_form_adm = ufl.inner(dt * b[i] * adj, v) * DX + sum(
            [dt * float(a[j, i]) * safe_action(safe_adjoint(derivative(
                forward_forms[j], y_)), kbar[j]) for j in range(i, size)])
        if explicit:
            stage_forms_adm = [stage_form_adm]
            jacobian_indices.append(-1)
        else:
            # Create a F=0 form and differentiate it
            stage_form_adm -= ufl.inner(kbar[i], v) * DX
            stage_forms_adm = [stage_form_adm, derivative(stage_form_adm, kbari)]
            jacobian_indices.append(1)

        ufl_stage_forms.append(stage_forms_adm)
        dolfin_stage_forms.append([Form(form) for form in stage_forms_adm])
        stage_solutions.append(kbari)

    # Only one last stage
    if len(b.shape) == 1:
        last_stage = Form(ufl.inner(adj, v) * DX + sum(
            [safe_action(safe_adjoint(derivative(forward_forms[i], y_)), kbar[i])
             for i in range(size)]))
    else:
        raise Exception("Not sure what to do here")

    human_form = "unimplemented"

    return ufl_stage_forms, dolfin_stage_forms, jacobian_indices, last_stage, \
        stage_solutions, dt, human_form, adj


class MultiStageScheme(cpp.multistage.MultiStageScheme):
    """Base class for all MultiStageSchemes

    """
    def __init__(self, rhs_form, ufl_stage_forms,
                 dolfin_stage_forms, last_stage, stage_solutions,
                 solution, time, dt, dt_stage_offsets, jacobian_indices, order,
                 name, human_form, bcs, contraction=None):

        # Store Python data
        self._rhs_form = rhs_form
        self._ufl_stage_forms = ufl_stage_forms
        self._dolfin_stage_forms = dolfin_stage_forms
        self._t = time
        self._bcs = bcs
        self._dt = dt
        self._last_stage = last_stage
        self._solution = solution
        self._stage_solutions = stage_solutions
        self._order = order
        self.jacobian_indices = jacobian_indices
        self.contraction = contraction

        # Pass args to C++ constructor
        stage_solutions = [s.cpp_object() for s in stage_solutions]

        cpp.multistage.MultiStageScheme.__init__(self,
                                                 dolfin_stage_forms,
                                                 last_stage,
                                                 stage_solutions,
                                                 solution.cpp_object(),
                                                 time.cpp_object(),
                                                 dt.cpp_object(),
                                                 dt_stage_offsets,
                                                 jacobian_indices,
                                                 order,
                                                 self.__class__.__name__,
                                                 human_form, bcs)

    def rhs_form(self):
        "Return the original rhs form"
        return self._rhs_form

    def ufl_stage_forms(self):
        "Return the ufl stage forms"
        return self._ufl_stage_forms

    def dolfin_stage_forms(self):
        "Return the dolfin stage forms"
        return self._dolfin_stage_forms

    def t(self):
        "Return the Constant used to describe time in the MultiStageScheme"
        return self._t

    def dt(self):
        "Return the Constant used to describe time in the MultiStageScheme"
        return self._dt

    def solution(self):
        "Return the solution Function"
        return self._solution

    def last_stage(self):
        "Return the form describing the last stage"
        return self._last_stage

    def stage_solutions(self):
        "Return the stage solutions"
        return self._stage_solutions

    def to_tlm(self, perturbation):
        raise NotImplementedError("'to_tlm:' implement in derived classes")

    def to_adm(self, perturbation):
        raise NotImplementedError("'to_adm:' implement in derived classes")


class ButcherMultiStageScheme(MultiStageScheme):
    """Base class for all MultiStageSchemes

    """
    def __init__(self, rhs_form, solution, time, bcs, a, b, c, order,
                 generator=_butcher_scheme_generator):
        bcs = bcs or []
        time = time or Constant(0.0)
        ufl_stage_forms, dolfin_stage_forms, jacobian_indices, last_stage, \
            stage_solutions, dt, human_form, contraction = \
            generator(a, b, c, time, solution, rhs_form)

        # Store data
        self.a = a
        self.b = b
        self.c = c

        MultiStageScheme.__init__(self, rhs_form, ufl_stage_forms,
                                  dolfin_stage_forms, last_stage,
                                  stage_solutions, solution, time, dt,
                                  c, jacobian_indices, order,
                                  self.__class__.__name__, human_form,
                                  bcs, contraction)

    def to_tlm(self, perturbation):
        r"""Return another MultiStageScheme that implements the tangent
        linearisation of the ODE solver.

        This takes \dot{y_n} (the derivative of y_n with respect to a
        parameter) and computes \dot{y_n+1} (the derivative of y_n+1
        with respect to that parameter).

        """

        generator = functools.partial(_butcher_scheme_generator_tlm,
                                      perturbation=perturbation)
        new_solution = self._solution.copy()
        new_form = ufl.replace(self._rhs_form, {self._solution: new_solution})
        return ButcherMultiStageScheme(new_form, new_solution,
                                       self._t, self._bcs, self.a,
                                       self.b, self.c, self._order,
                                       generator=generator)

    def to_adm(self, adj):
        r"""Return another MultiStageScheme that implements the adjoint
        linearisation of the ODE solver.

        This takes \bar{y_n+1} (the derivative of a functional J with
        respect to y_n+1) and computes \bar{y_n} (the derivative of J
        with respect to y_n).

        """

        generator = functools.partial(_butcher_scheme_generator_adm, adj=adj)
        new_solution = self._solution.copy()
        new_form = ufl.replace(self._rhs_form, {self._solution: new_solution})
        return ButcherMultiStageScheme(new_form, new_solution, self._t, self._bcs,
                                       self.a, self.b, self.c, self._order,
                                       generator=generator)


class ERK1(ButcherMultiStageScheme):
    """Explicit first order Scheme"""
    def __init__(self, rhs_form, solution, t=None, bcs=None):
        a = np.array([0.])
        b = np.array([1.])
        c = np.array([0.])
        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs,
                                         a, b, c, 1)


class BDF1(ButcherMultiStageScheme):
    """Implicit first order scheme"""
    def __init__(self, rhs_form, solution, t=None, bcs=None):
        a = np.array([1.])
        b = np.array([1.])
        c = np.array([1.])
        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs,
                                         a, b, c, 1)


class ExplicitMidPoint(ButcherMultiStageScheme):
    """Explicit 2nd order scheme"""
    def __init__(self, rhs_form, solution, t=None, bcs=None):

        a = np.array([[0, 0], [0.5, 0.0]])
        b = np.array([0., 1])
        c = np.array([0, 0.5])
        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs,
                                         a, b, c, 2)


class CN2(ButcherMultiStageScheme):
    """Semi-implicit 2nd order scheme"""
    def __init__(self, rhs_form, solution, t=None, bcs=None):
        a = np.array([[0, 0], [0.5, 0.5]])
        b = np.array([0.5, 0.5])
        c = np.array([0, 1.0])

        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs,
                                         a, b, c, 2)


class ERK4(ButcherMultiStageScheme):
    """Explicit 4th order scheme"""
    def __init__(self, rhs_form, solution, t=None, bcs=None):
        a = np.array([[0, 0, 0, 0],
                      [0.5, 0, 0, 0],
                      [0, 0.5, 0, 0],
                      [0, 0, 1, 0]])
        b = np.array([1. / 6, 1. / 3, 1. / 3, 1. / 6])
        c = np.array([0, 0.5, 0.5, 1])
        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs,
                                         a, b, c, 4)


class ESDIRK3(ButcherMultiStageScheme):
    """Explicit implicit 3rd order scheme

    See also "Singly diagonally implicit Runge–Kutta methods with an
    explicit first stage" by A Kværnø - BIT Numerical Mathematics,
    2004 (p.497)

    """
    def __init__(self, rhs_form, solution, t=None, bcs=None):
        a = np.array([[0.000000000000000, 0.000000000000000, 0.000000000000000, 0.00000000000000],
                      [0.435866521500000, 0.435866521500000, 0.000000000000000, 0.00000000000000],
                      [0.490563388419108, 0.073570090080892, 0.435866521500000, 0.00000000000000],
                      [0.308809969973036, 1.490563388254108, -1.235239879727145, 0.435866521500000]])
        b = a[-1, :].copy()
        c = a.sum(1)
        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs, a, b, c, 3)


class ESDIRK4(ButcherMultiStageScheme):
    """Explicit implicit 4rd order scheme

    See also "Singly diagonally implicit Runge–Kutta methods with an
    explicit first stage" by A Kværnø - BIT Numerical Mathematics,
    2004 (p.498)

    """
    def __init__(self, rhs_form, solution, t=None, bcs=None):
        a = np.array([[0.000000000000000, 0.0000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000],
                      [0.435866521500000, 0.4358665215000000, 0.000000000000000, 0.000000000000000, 0.000000000000000],
                      [0.140737774731968, -0.108365551378832, 0.435866521500000, 0.000000000000000, 0.000000000000000],
                      [0.102399400616089, -0.376878452267324, 0.838612530151233, 0.435866521500000, 0.000000000000000],
                      [0.157024897860995, 0.1173304413577680, 0.616678030391680, -0.326899891110444, 0.435866521500000]])

        b = a[-1, :].copy()
        c = a.sum(1)
        ButcherMultiStageScheme.__init__(self, rhs_form, solution, t, bcs, a, b, c, 4)


# Aliases
CrankNicolson = CN2
ExplicitEuler = ERK1
ForwardEuler = ERK1
ImplicitEuler = BDF1
BackwardEuler = BDF1
ERK = ERK1
RK4 = ERK4

__all__ = [name for name, attr in list(globals().items())
           if isinstance(attr, type) and issubclass(attr, MultiStageScheme)]

__all__.append("MultiStageScheme")
