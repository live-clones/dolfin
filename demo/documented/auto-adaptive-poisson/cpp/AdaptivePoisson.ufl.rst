UFL input for the auto adaptive Poisson problem
===============================================

UFL code::

  from ufl_legacy import (Coefficient, dot, ds, dx, FiniteElement, grad,
                   TestFunction, TrialFunction, triangle)

  element = FiniteElement("CG", triangle, 1)
  u = TrialFunction(element)
  v = TestFunction(element)

  f = Coefficient(element)
  g = Coefficient(element)

  a = dot(grad(u), grad(v))*dx()
  L = f*v*dx() + g*v*ds()
  M = u*dx()

Before the form file can be used in the C++ program, it must be
compiled using FFC by running (on the command-line):

.. code-block:: sh

   ffc -l dolfin -e AdaptivePoisson.ufl

Parameter ``-e`` ensures that a code for forms used for error control
is generated.
