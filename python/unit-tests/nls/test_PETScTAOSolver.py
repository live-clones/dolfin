#!/usr/bin/env py.test

"""Unit test for the PETSc TAO solver"""

# Copyright (C) 2014 Tianyi Li
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
# First added:  2014-07-19
# Last changed: 2014-10-14

from __future__ import print_function

"""This demo uses PETSc TAO solver for nonlinear (bound-constrained)
optimisation problems to solve the same minimization problem for testing
TAOLinearBoundSolver."""

from dolfin import *
import pytest

from dolfin_utils.test import *

backend = set_parameters_fixture("linear_algebra_backend", ["PETSc"])

@pytest.mark.xfail
@skip_if_not_PETSc
def test_tao_linear_bound_solver(backend):
    "Test PETScTAOSolver"

    # Create mesh and define function space
    Lx, Ly = 1.0, 0.1
    mesh = RectangleMesh(Point(0, 0), Point(Lx, Ly), 100, 10)
    V = FunctionSpace(mesh, "Lagrange", 1)

    # Define Dirichlet boundaries
    def left(x,on_boundary):
        return on_boundary and near(x[0], 0.0)

    def right(x,on_boundary):
        return on_boundary and near(x[0], Lx)

    # Define boundary conditions
    zero = Constant(0.0)
    one  = Constant(1.0)
    bc_l = DirichletBC(V, zero, left)
    bc_r = DirichletBC(V, one, right)
    bcs = [bc_l, bc_r]

    # Define the variational problem
    u = Function(V)
    du = TrialFunction(V)
    v = TestFunction(V)
    cv = Constant(3.0/4.0)
    ell = Constant(0.5)  # This should be smaller than Lx
    energy = cv*(ell/2.0*inner(grad(u), grad(u))*dx + 2.0/ell*u*dx)
    grad_energy = derivative(energy, u, v)
    H_energy = derivative(grad_energy, u, du)

    # Define the lower and upper bounds
    lb = interpolate(Constant(0.0), V)
    ub = interpolate(Constant(1.0), V)

    # Apply BC to the lower and upper bounds
    for bc in bcs:
        bc.apply(lb.vector())
        bc.apply(ub.vector())

    # Define the minimisation problem by using OptimisationProblem class
    class TestProblem(OptimisationProblem):

        def __init__(self, u, energy, grad_energy, H_energy):
            OptimisationProblem.__init__(self)
            self.u = u
            self.energy = energy
            self.grad_energy = grad_energy
            self.H_energy = H_energy

        # Objective function
        def f(self, x):
            print("WWWWWWWWW")
            self.u.vector()[:] = x
            return assemble(self.energy)

        # Gradient of the objective function
        def F(self, b, x):
            self.u.vector()[:] = x
            assemble(self.grad_energy, tensor=b)

        # Hessian of the objective function
        def J(self, A, x):
            self.u.vector()[:] = x
            assemble(self.H_energy, tensor=A)

    # Create the PETScTAOSolver
    solver = PETScTAOSolver()

    # Set some parameters
    solver.parameters["monitor_convergence"] = True
    solver.parameters["report"] = True

    p = TestProblem(u, energy, grad_energy, H_energy)
    solver.solve(p, u.vector(), lb.vector(), ub.vector())
    #solver.solve(TestProblem(), u.vector(), lb.vector(), ub.vector())

    # Verify that energy(u) = Ly
    #assert round(assemble(energy) - Ly, 4) == 0
