# -*- coding: utf-8 -*-
"""This module defines different MultiStageScheme classes based on the
Rush Larsen expicit integration schemes which can be passed to a
PointIntegralSolver

"""

# Copyright (C) 2014 Johan Hake
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

import functools
import ufl_legacy as ufl

from dolfin.function.constant import Constant
from dolfin.function.function import Function
from dolfin.fem.formmanipulations import derivative
from dolfin.multistage.factorize import extract_tested_expressions
from dolfin.fem.form import Form
from dolfin.multistage.multistagescheme import (MultiStageScheme,
                                                _check_form, _time_dependent_expressions,
                                                _replace_dict_time_dependent_expression, safe_action,
                                                safe_adjoint)
from dolfin import DOLFIN_EPS
from dolfin import TrialFunction

import ufl_legacy.algorithms
from ufl_legacy.algorithms import (expand_derivatives, expand_indices,
                                   extract_coefficients)


def _rush_larsen_step(rhs_exprs, diff_rhs_exprs, linear_terms,
                      system_size, y0, stage_solution, dt, time, a, c,
                      v, DX, time_dep_expressions):

    # If we need to replace the original solution with stage solution
    repl = None
    if stage_solution is not None:

        if system_size > 1:
            repl = {y0: stage_solution}
        else:
            repl = {y0[0]: stage_solution}

    # If we have time dependent expressions
    if time_dep_expressions and abs(float(c)) > DOLFIN_EPS:
        time_ = time
        time = time + dt * float(c)

        repl.update(_replace_dict_time_dependent_expression(time_dep_expressions,
                                                            time_, dt, float(c)))
        repl[time_] = time

    # If all terms are linear (using generalized=True) we add a safe
    # guard to the linearized term. See below
    safe_guard = sum(linear_terms) == system_size

    # Add componentwise contribution to rl form
    rl_ufl_form = ufl.zero()
    num_rl_steps = 0
    for ind in range(system_size):

        # forward euler step
        fe_du_i = rhs_exprs[ind] * dt * float(a)

        # If exact integration
        if linear_terms[ind]:
            num_rl_steps += 1

            # Rush Larsen step Safeguard the divisor: let's hope
            # diff_rhs_exprs[ind] is never 1.0e-16!  Let's get rid of
            # this when the conditional fixes land properly in UFL.
            eps = Constant(1.0e-16)
            rl_du_i = rhs_exprs[ind] / (diff_rhs_exprs[ind] + eps) * (
                ufl.exp(diff_rhs_exprs[ind] * dt) - 1.0)

            # If safe guard
            if safe_guard:
                du_i = ufl.conditional(ufl.lt(abs(diff_rhs_exprs[ind]), 1e-8), fe_du_i, rl_du_i)
            else:
                du_i = rl_du_i
        else:
            du_i = fe_du_i

        # If we should replace solution in form with stage solution
        if repl:
            du_i = ufl.replace(du_i, repl)

        rl_ufl_form += (y0[ind] + du_i) * v[ind]

    return rl_ufl_form * DX


def _find_linear_terms(rhs_exprs, u):
    """Help function that takes a list of rhs expressions and return a
    list of bools determining what component, rhs_exprs[i], is linear
    wrt u[i].

    """

    uu = [Constant(1.0) for _ in rhs_exprs]
    if len(rhs_exprs) > 1:
        repl = {u: ufl.as_vector(uu)}
    else:
        repl = {u: uu[0]}

    linear_terms = []
    for i, ui in enumerate(uu):
        comp_i_s = expand_indices(ufl.replace(rhs_exprs[i], repl))
        linear_terms.append(ui in extract_coefficients(comp_i_s) and
                            ui not in extract_coefficients(
                                expand_derivatives(ufl.diff(comp_i_s, ui))))
    return linear_terms


def _rush_larsen_scheme_generator(rhs_form, solution, time, order, generalized):
    """Generates a list of forms and solutions for a given Butcher
    tableau

    *Arguments*
        rhs_form (ufl.Form)
            A UFL form representing the rhs for a time differentiated equation
        solution (_Function_)
            The prognostic variable
        time (_Constant_)
            A Constant holding the time at the start of the time step
        order (int)
            The order of the scheme
        generalized (bool)
            If True generate a generalized Rush Larsen scheme, linearizing all
            components.

    """

    DX = _check_form(rhs_form)

    if DX != ufl.dP:
        raise TypeError("Expected a form with a Pointintegral.")

    # Create time step
    dt = Constant(0.1)

    # Get test function
    #    arguments = rhs_form.arguments()
    #    coefficients = rhs_form.coefficients()

    # Get time dependent expressions
    time_dep_expressions = _time_dependent_expressions(rhs_form, time)

    # Extract rhs expressions from form
    rhs_integrand = rhs_form.integrals()[0].integrand()
    rhs_exprs, v = extract_tested_expressions(rhs_integrand)
    vector_rhs = len(v.ufl_shape) > 0 and v.ufl_shape[0] > 1

    system_size = v.ufl_shape[0] if vector_rhs else 1

    # Fix for indexing of v for scalar expressions
    v = v if vector_rhs else [v]

    # Extract linear terms if not using generalized Rush Larsen
    if not generalized:
        linear_terms = _find_linear_terms(rhs_exprs, solution)
    else:
        linear_terms = [True for _ in range(system_size)]

    # Wrap the rhs expressions into a ufl vector type
    rhs_exprs = ufl.as_vector([rhs_exprs[i] for i in range(system_size)])
    rhs_jac = ufl.diff(rhs_exprs, solution)

    # Takes time!
    if vector_rhs:
        diff_rhs_exprs = [expand_indices(expand_derivatives(rhs_jac[ind, ind]))
                          for ind in range(system_size)]
    else:
        diff_rhs_exprs = [expand_indices(expand_derivatives(rhs_jac[0]))]
        solution = [solution]

    ufl_stage_forms = []
    dolfin_stage_forms = []
    dt_stage_offsets = []

    # Stage solutions (3 per order rhs, linearized, and final step)
    # If 2nd order the final step for 1 step is a stage
    if order == 1:
        stage_solutions = []
        rl_ufl_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs, linear_terms,
                                        system_size, solution, None, dt, time, 1.0,
                                        0.0, v, DX, time_dep_expressions)
    elif order == 2:

        # Stage solution for order 2
        if vector_rhs:
            stage_solutions = [Function(solution.function_space(), name="y_1/2")]
        else:
            stage_solutions = [Function(solution[0].function_space(), name="y_1/2")]

        stage_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                       linear_terms, system_size,
                                       solution, None, dt, time, 0.5,
                                       0.0, v, DX,
                                       time_dep_expressions)

        rl_ufl_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                        linear_terms, system_size,
                                        solution, stage_solutions[0],
                                        dt, time, 1.0, 0.5, v, DX,
                                        time_dep_expressions)

        ufl_stage_forms.append([stage_form])
        dolfin_stage_forms.append([Form(stage_form)])

    # Get last stage form
    last_stage = Form(rl_ufl_form)

    human_form = "%srush larsen %s" % ("generalized " if generalized else "",
                                       str(order))

    return rhs_form, linear_terms, ufl_stage_forms, dolfin_stage_forms, last_stage, \
        stage_solutions, dt, dt_stage_offsets, human_form, None


def _rush_larsen_scheme_generator_tlm(rhs_form, solution, time, order,
                                      generalized, perturbation):
    """Generates a list of forms and solutions for the tangent
    linearisation of the Rush-Larsen scheme

    *Arguments*
        rhs_form (ufl.Form)
            A UFL form representing the rhs for a time differentiated equation
        solution (_Function_)
            The prognostic variable
        time (_Constant_)
            A Constant holding the time at the start of the time step
        order (int)
            The order of the scheme
        generalized (bool)
            If True generate a generalized Rush Larsen scheme, linearizing all
            components.
        perturbation (Function)
            The vector on which we compute the tangent linear action.

    """

    DX = _check_form(rhs_form)

    if DX != ufl.dP:
        raise TypeError("Expected a form with a Pointintegral.")

    # Create time step
    dt = Constant(0.1)

    # Get test function
    #    arguments = rhs_form.arguments()
    #    coefficients = rhs_form.coefficients()

    # Get time dependent expressions
    time_dep_expressions = _time_dependent_expressions(rhs_form, time)

    # Extract rhs expressions from form
    rhs_integrand = rhs_form.integrals()[0].integrand()
    rhs_exprs, v = extract_tested_expressions(rhs_integrand)
    vector_rhs = len(v.ufl_shape) > 0 and v.ufl_shape[0] > 1

    system_size = v.ufl_shape[0] if vector_rhs else 1

    # Fix for indexing of v for scalar expressions
    v = v if vector_rhs else [v]

    # Extract linear terms if not using generalized Rush Larsen
    if not generalized:
        linear_terms = _find_linear_terms(rhs_exprs, solution)
    else:
        linear_terms = [True for _ in range(system_size)]

    # Wrap the rhs expressions into a ufl vector type
    rhs_exprs = ufl.as_vector([rhs_exprs[i] for i in range(system_size)])
    rhs_jac = ufl.diff(rhs_exprs, solution)

    # Takes time!
    if vector_rhs:
        diff_rhs_exprs = [expand_indices(expand_derivatives(rhs_jac[ind, ind]))
                          for ind in range(system_size)]
        soln = solution
    else:
        diff_rhs_exprs = [expand_indices(expand_derivatives(rhs_jac[0]))]
        solution = [solution]
        soln = solution[0]

    ufl_stage_forms = []
    dolfin_stage_forms = []
    dt_stage_offsets = []
    trial = TrialFunction(soln.function_space())

    # Stage solutions (3 per order rhs, linearized, and final step)
    # If 2nd order the final step for 1 step is a stage
    if order == 1:
        stage_solutions = []

        # Fetch the original step
        rl_ufl_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                        linear_terms, system_size,
                                        solution, None, dt, time, 1.0,
                                        0.0, v, DX,
                                        time_dep_expressions)

        rl_ufl_form = safe_action(derivative(rl_ufl_form, soln, trial), perturbation)

    elif order == 2:
        # Stage solution for order 2
        fn_space = soln.function_space()

        stage_solutions = [Function(fn_space, name="y_1/2"),
                           Function(fn_space, name="y_dot_1/2"),
                           Function(fn_space, name="y_1")]

        y_half_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                        linear_terms, system_size,
                                        solution, None, dt, time, 0.5,
                                        0.0, v, DX,
                                        time_dep_expressions)

        y_dot_half_form = safe_action(derivative(y_half_form, soln, trial), perturbation)

        y_one_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                       linear_terms, system_size,
                                       solution, stage_solutions[0],
                                       dt, time, 1.0, 0.5, v, DX,
                                       time_dep_expressions)

        rl_ufl_form = safe_action(derivative(y_one_form, soln, trial), perturbation) + \
            safe_action(derivative(y_one_form, stage_solutions[0], trial), stage_solutions[1])

        ufl_stage_forms.append([y_half_form])
        ufl_stage_forms.append([y_dot_half_form])
        ufl_stage_forms.append([y_one_form])
        dolfin_stage_forms.append([Form(y_half_form)])
        dolfin_stage_forms.append([Form(y_dot_half_form)])
        dolfin_stage_forms.append([Form(y_one_form)])

    # Get last stage form
    last_stage = Form(rl_ufl_form)

    human_form = "%srush larsen %s" % ("generalized " if generalized else "",
                                       str(order))

    return rhs_form, linear_terms, ufl_stage_forms, dolfin_stage_forms, last_stage, \
        stage_solutions, dt, dt_stage_offsets, human_form, perturbation


def _rush_larsen_scheme_generator_adm(rhs_form, solution, time, order,
                                      generalized, perturbation):
    """Generates a list of forms and solutions for the adjoint
    linearisation of the Rush-Larsen scheme

    *Arguments*
        rhs_form (ufl.Form)
            A UFL form representing the rhs for a time differentiated equation
        solution (_Function_)
            The prognostic variable
        time (_Constant_)
            A Constant holding the time at the start of the time step
        order (int)
            The order of the scheme
        generalized (bool)
            If True generate a generalized Rush Larsen scheme, linearizing all
            components.
        perturbation (Function)
            The vector on which we compute the adjoint action.

    """

    DX = _check_form(rhs_form)

    if DX != ufl.dP:
        raise TypeError("Expected a form with a Pointintegral.")

    # Create time step
    dt = Constant(0.1)

    # Get test function
    #    arguments = rhs_form.arguments()
    #    coefficients = rhs_form.coefficients()

    # Get time dependent expressions
    time_dep_expressions = _time_dependent_expressions(rhs_form, time)

    # Extract rhs expressions from form
    rhs_integrand = rhs_form.integrals()[0].integrand()
    rhs_exprs, v = extract_tested_expressions(rhs_integrand)
    vector_rhs = len(v.ufl_shape) > 0 and v.ufl_shape[0] > 1

    system_size = v.ufl_shape[0] if vector_rhs else 1

    # Fix for indexing of v for scalar expressions
    v = v if vector_rhs else [v]

    # Extract linear terms if not using generalized Rush Larsen
    if not generalized:
        linear_terms = _find_linear_terms(rhs_exprs, solution)
    else:
        linear_terms = [True for _ in range(system_size)]

    # Wrap the rhs expressions into a ufl vector type
    rhs_exprs = ufl.as_vector([rhs_exprs[i] for i in range(system_size)])
    rhs_jac = ufl.diff(rhs_exprs, solution)

    # Takes time!
    if vector_rhs:
        diff_rhs_exprs = [expand_indices(expand_derivatives(rhs_jac[ind, ind]))
                          for ind in range(system_size)]
        soln = solution
    else:
        diff_rhs_exprs = [expand_indices(expand_derivatives(rhs_jac[0]))]
        solution = [solution]
        soln = solution[0]

    ufl_stage_forms = []
    dolfin_stage_forms = []
    dt_stage_offsets = []
    trial = TrialFunction(soln.function_space())

    # Stage solutions (3 per order rhs, linearized, and final step)
    # If 2nd order the final step for 1 step is a stage
    if order == 1:
        stage_solutions = []

        # Fetch the original step
        rl_ufl_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                        linear_terms, system_size,
                                        solution, None, dt, time, 1.0,
                                        0.0, v, DX,
                                        time_dep_expressions)

        # If this is commented out, we don't get NaNs.  Yhy is
        # solution a list of length zero anyway?
        rl_ufl_form = safe_action(safe_adjoint(derivative(rl_ufl_form, soln, trial)),
                                  perturbation)

    elif order == 2:
        # Stage solution for order 2
        fn_space = soln.function_space()

        stage_solutions = [Function(fn_space, name="y_1/2"),
                           Function(fn_space, name="y_1"),
                           Function(fn_space, name="y_bar_1/2")]

        y_half_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                        linear_terms, system_size,
                                        solution, None, dt, time, 0.5,
                                        0.0, v, DX,
                                        time_dep_expressions)

        y_one_form = _rush_larsen_step(rhs_exprs, diff_rhs_exprs,
                                       linear_terms, system_size,
                                       solution, stage_solutions[0],
                                       dt, time, 1.0, 0.5, v, DX,
                                       time_dep_expressions)

        y_bar_half_form = safe_action(safe_adjoint(derivative(y_one_form,
                                                              stage_solutions[0], trial)), perturbation)

        rl_ufl_form = safe_action(safe_adjoint(derivative(y_one_form, soln, trial)), perturbation) + \
            safe_action(safe_adjoint(derivative(y_half_form, soln, trial)), stage_solutions[2])

        ufl_stage_forms.append([y_half_form])
        ufl_stage_forms.append([y_one_form])
        ufl_stage_forms.append([y_bar_half_form])
        dolfin_stage_forms.append([Form(y_half_form)])
        dolfin_stage_forms.append([Form(y_one_form)])
        dolfin_stage_forms.append([Form(y_bar_half_form)])

    # Get last stage form
    last_stage = Form(rl_ufl_form)

    human_form = "%srush larsen %s" % ("generalized " if generalized else "",
                                       str(order))

    return rhs_form, linear_terms, ufl_stage_forms, dolfin_stage_forms, last_stage, \
        stage_solutions, dt, dt_stage_offsets, human_form, perturbation


class RushLarsenScheme(MultiStageScheme):
    def __init__(self, rhs_form, solution, time, order, generalized,
                 generator=_rush_larsen_scheme_generator):

        self._rhs_form = rhs_form
        self._solution = solution
        self._t = time
        self._order = order
        self._generalized = generalized

        # FIXME: What with bcs?
        bcs = []
        time = time or Constant(0.0)
        if order not in [1, 2]:
            raise ValueError("Expected order to be either 1 or 2")

        rhs_form, ufl_stage_forms, linear_terms, dofin_stage_forms, last_stage, \
            stage_solutions, dt, dt_stage_offsets, human_form, contraction = \
            generator(rhs_form, solution, time, order, generalized)

        self.linear_terms = linear_terms

        # All stages are explicit
        jacobian_indices = [-1 for _ in stage_solutions]

        # Highjack a and b to hold parameters for RushLarsen scheme generation
        MultiStageScheme.__init__(self, rhs_form, ufl_stage_forms,
                                  dofin_stage_forms, last_stage,
                                  stage_solutions, solution, time, dt,
                                  dt_stage_offsets, jacobian_indices,
                                  order, self.__class__.__name__,
                                  human_form, bcs, contraction)

    def to_tlm(self, perturbation):
        r"""Return another RushLarsenScheme that implements the tangent
        linearisation of the ODE solver.

        This takes \dot{y_n} (the derivative of y_n with respect to a
        parameter) and computes \dot{y_n+1} (the derivative of y_n+1
        with respect to that parameter).

        """

        generator = functools.partial(_rush_larsen_scheme_generator_tlm,
                                      perturbation=perturbation)
        new_solution = self._solution.copy()
        new_form = ufl.replace(self._rhs_form, {self._solution: new_solution})
        return RushLarsenScheme(new_form, new_solution, self._t, self._order,
                                self._generalized, generator=generator)

    def to_adm(self, perturbation):
        r"""Return another RushLarsenScheme that implements the adjoint
        linearisation of the ODE solver.

        This takes \bar{y_n+1} (the derivative of y_n+1 with respect to a
        parameter) and computes \bar{y_n} (the derivative of y_n
        with respect to that parameter).

        """

        generator = functools.partial(_rush_larsen_scheme_generator_adm,
                                      perturbation=perturbation)
        new_solution = self._solution.copy()
        new_form = ufl.replace(self._rhs_form, {self._solution: new_solution})
        return RushLarsenScheme(new_form, new_solution, self._t,
                                self._order, self._generalized,
                                generator=generator)


class RL1(RushLarsenScheme):
    """First order Rush Larsen Scheme

    """
    def __init__(self, rhs_form, solution, t=None):
        RushLarsenScheme.__init__(self, rhs_form, solution, t, 1, False)


class RL2(RushLarsenScheme):
    """Second order Rush Larsen Scheme

    """
    def __init__(self, rhs_form, solution, t=None):
        RushLarsenScheme.__init__(self, rhs_form, solution, t, 2, generalized=False)


class GRL1(RushLarsenScheme):
    """
    First order generalized Rush Larsen Scheme
    """
    def __init__(self, rhs_form, solution, t=None):
        RushLarsenScheme.__init__(self, rhs_form, solution, t, 1, True)


class GRL2(RushLarsenScheme):
    """Second order generalized Rush Larsen Scheme

    """
    def __init__(self, rhs_form, solution, t=None):
        RushLarsenScheme.__init__(self, rhs_form, solution, t, 2, True)


__all__ = [name for name, attr in list(globals().items())
           if isinstance(attr, type) and issubclass(attr, MultiStageScheme)]

__all__.append("MultiStageScheme")
