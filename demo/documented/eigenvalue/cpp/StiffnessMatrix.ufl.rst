UFL input for the Poisson bilinear form
=======================================

The bilinear form for a stiffness matrix (Poisson)::

  from ufl import (dot, dx, FiniteElement, grad,
                   TestFunction, TrialFunction, tetrahedron)

  element = FiniteElement("Lagrange", tetrahedron, 1)

  v = TestFunction(element)
  u = TrialFunction(element)

  a = dot(grad(v), grad(u))*dx
