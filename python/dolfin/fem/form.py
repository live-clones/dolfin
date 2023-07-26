# -*- coding: utf-8 -*-
"""FIXME: Add description"""

# Copyright (C) 2017 Chris N. Richardson and Garth N. Wells
#
# Distributed under the terms of the GNU Lesser Public License (LGPL),
# either version 3 of the License, or (at your option) any later
# version.

import ufl_legacy as ufl
import dolfin.cpp as cpp
from dolfin.jit.jit import dolfin_pc, ffc_jit


class Form(cpp.fem.Form):
    def __init__(self, form, **kwargs):

        # Check form argument
        if not isinstance(form, ufl.Form):
            raise RuntimeError("Expected a ufl.Form.")

        sd = form.subdomain_data()
        self.subdomains, = list(sd.values())  # Assuming single domain
        domain, = list(sd.keys())  # Assuming single domain
        mesh = domain.ufl_cargo()
        if isinstance(mesh, cpp.mesh.MultiMesh):
            mesh = mesh.part(0)

        # Having a mesh in the form is a requirement
        if mesh is None:
            raise RuntimeError("Expecting to find a Mesh in the form.")

        form_compiler_parameters = kwargs.pop("form_compiler_parameters", None)

        # Add DOLFIN include paths (just the Boost path for special
        # math functions is really required)
        # FIXME: move getting include paths to elsewhere
        # FIXME: How to add these path is form_compiler parameters is input,
        # and a dolfin::Parameters
        if form_compiler_parameters is None:
            form_compiler_parameters = {"external_include_dirs": dolfin_pc["include_dirs"]}

        ufc_form = ffc_jit(form, form_compiler_parameters=form_compiler_parameters,
                           mpi_comm=mesh.mpi_comm())
        ufc_form = cpp.fem.make_ufc_form(ufc_form[0])

        # TO BE CHECKED
        self.function_spaces = kwargs.get("function_spaces")

        # Extraction of functionspaces contained in a MultiMeshFunctionSpace
        if not self.function_spaces:
            self.function_spaces = [func.ufl_function_space()._cpp_object for func in form.arguments()]

        cpp.fem.Form.__init__(self, ufc_form, self.function_spaces)

        original_coefficients = form.coefficients()
        self.coefficients = []
        for i in range(self.num_coefficients()):
            j = self.original_coefficient_position(i)
            self.coefficients.append(original_coefficients[j]._cpp_object)

        # Type checking coefficients
        if not all(isinstance(c, (cpp.function.GenericFunction, cpp.function.MultiMeshFunction))
                   for c in self.coefficients):
            coefficient_error = "Error while extracting coefficients. "
            raise TypeError(coefficient_error +
                            "Either provide a dict of cpp.function.GenericFunctions, " +
                            "or use Function to define your form.")

        for i in range(self.num_coefficients()):
            if isinstance(self.coefficients[i], cpp.function.GenericFunction):
                self.set_coefficient(i, self.coefficients[i])

        # Attach mesh :
        # - because function spaces and coefficients may be empty lists
        # - because function spaces can be built from different meshes
        self.set_mesh(mesh)

        # Attach subdomains to C++ Form if we have them
        subdomains = self.subdomains.get("cell")
        if subdomains is not None:
            self.set_cell_domains(subdomains)
        subdomains = self.subdomains.get("exterior_facet")
        if subdomains is not None:
            self.set_exterior_facet_domains(subdomains)
        subdomains = self.subdomains.get("interior_facet")
        if subdomains is not None:
            self.set_interior_facet_domains(subdomains)
        subdomains = self.subdomains.get("vertex")
        if subdomains is not None:
            self.set_vertex_domains(subdomains)
