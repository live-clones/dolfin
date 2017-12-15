# -*- coding: utf-8 -*-

# Copyright (C) 2017 Jørgen S. Dokken
#
# Distributed under the terms of the GNU Lesser Public License (LGPL),
# either version 3 of the License, or (at your option) any later
# version.

import ufl
import dolfin.cpp as cpp
from dolfin.function.functionspace import FunctionSpace
from six import string_types

class MultiMeshFunctionSpace(cpp.function.MultiMeshFunctionSpace):
    def __init__(self, *args, **kwargs):
        """Create multimesh finite element function space.
        
        *Arguments*
        multimesh
        a :py:class:`MultiMesh <dolfin.cpp.mesh.MultiMesh>`.
        family
        a string specifying the element family,
            see :py:class:`FunctionSpace
        <dolfin.functions.functionspace.FunctionSpace>`
        for alternatives.
        
            This argument may also be a `FiniteElement`, in
        which case the `degree` argument should not be
        specified.
        degree
        the degree of the element.
        
        *Example of usage*
        
        .. code-block:: python
        
        V = MultiMeshFunctionSpace(mesh, "CG", 1)
        
        element = FiniteElement("Lagrange", triangle, 1)
        V = MultiMeshFunctionSpace(mesh, element)
        """
        if len(args)==2:
            self.__init_from_ufl(*args, **kwargs)
        elif len(args)==3:
            self._init_convenience(*args, **kwargs)
        else:
            raise NotImplementedError

    def __init_from_ufl(self, multimesh, element):
        self.info = [element]
        if not isinstance(element, ufl.FiniteElementBase):
            cpp.dolfin_error("multimeshfunctionspace.py",
                             "create function space",
                             "Illegal argument, not a finite element: "
                             + str(element))

        # Create and add individual function spaces
        V = cpp.function.MultiMeshFunctionSpace(multimesh)
        V_parts = []
        for part in range(multimesh.num_parts()):
            V_part = FunctionSpace(multimesh.part(part), element)
            V_parts.append(V_part)
            V.add(V_part)

        # Build multimesh function space
        V.build()

        # Store full function spaces
        self._parts = V_parts
        self._cpp_object = V

    def _init_convenience(self, multimesh, family, degree):
        # Check arguments
        self.info = [family, degree]
        if not isinstance(family, string_types):
            cpp.dolfin_error("multimeshfunctionspace.py",
                             "create function space",
                             "Illegal argument for finite element family, not a string: " + str(family))
        if not isinstance(degree, int):
            cpp.dolfin_error("multimeshfunctionspace.py",
                             "create function space",
                             "Illegal argument for degree, not an integer: "
                             + str(degree))
        if not isinstance(multimesh, cpp.mesh.MultiMesh):
            cpp.dolfin_error("functionspace.py",
                             "create multimesh function space",
                             "Illegal argument, not a multimesh: " + str(multimesh))

        # Create UFL element
        mesh = multimesh.part(0)
        element = ufl.FiniteElement(family, mesh.ufl_cell(), degree)

        # Create and add individual function spaces
        V = cpp.function.MultiMeshFunctionSpace(multimesh)
        V_parts = []
        for part in range(multimesh.num_parts()):
            V_part = FunctionSpace(multimesh.part(part), element)
            V_parts.append(V_part)
            V.add(V_part)

        # Build multimesh function space
        V.build()
        
        # Store full function spaces
        self._parts = V_parts
        self._cpp_object = V
