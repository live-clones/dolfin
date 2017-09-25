# -*- coding: utf-8 -*-
"""FIXME: Add description"""

# Copyright (C) 2017 Chris N. Richardson and Garth N. Wells
#
# Distributed under the terms of the GNU Lesser Public License (LGPL),
# either version 3 of the License, or (at your option) any later
# version.

import pkgconfig

import ufl
import ffc

import sys
import traceback

import dolfin.cpp as cpp
from dolfin.parameter import ffc_default_parameters


class Form(cpp.fem.Form):
    def __init__(self, form, **kwargs):

        # Check form argument
        if not isinstance(form, ufl.Form):
            raise RuntimeError("Expected a ufl.Form.")

        sd = form.subdomain_data()
        self.subdomains, = list(sd.values())  # Assuming single domain
        domain, = list(sd.keys())  # Assuming single domain
        mesh = domain.ufl_cargo()

        # Having a mesh in the form is a requirement
        if mesh is None:
            raise RuntimeError("Expecting to find a Mesh in the form.")

        # FIXME: Have explicit kwarg in function signature
        form_compiler_parameters = kwargs.pop("form_compiler_parameters", None)

        # Prepare form compiler parameters with overrides from dolfin and kwargs
        p = ffc_default_parameters()
        p.update(cpp.parameter.parameters["form_compiler"])
        p.update(form_compiler_parameters or {})

        # Add DOLFIN include paths (just the Boost path for special
        # math functions is really required)
        # FIXME: move getting include paths to elsewhere
        # FIXME: do we really want to store include dirs as repr(str)?
        #        we would need to do
        #          p["external_include_dirs"] = repr(dirs + eval(p["external_include_dirs"]))
        #        which is a security risk
        dirs = repr(pkgconfig.parse('dolfin')["include_dirs"])
        assert p["external_include_dirs"] == "", "Don't know how to concat include dirs"
        p["external_include_dirs"] += dirs

        # FIXME: Cast Parameters to dict to avoid later problem in FFC; get rid of this block
        p_ = {}
        for k in p:
            p_[k] = p[k]
        p = p_

        # Execute!
        try:
            ufc_form = ffc.jit(form, p)
        except Exception as e:
            #tb_text = ''.join(traceback.format_exception(*sys.exc_info()))
            #raise RuntimeError("ffc.jit failed with message:\n%s" % (tb_text,))
            raise RuntimeError("ffc.jit failed; see message above")

        # Unpack result (ffc.jit returns different tuples based on input type)
        ufc_form = cpp.fem.make_ufc_form(ufc_form[0])

        function_spaces = [func.function_space()._cpp_object for func in form.arguments()]

        cpp.fem.Form.__init__(self, ufc_form, function_spaces)

        original_coefficients = form.coefficients()
        self.coefficients = []
        for i in range(self.num_coefficients()):
            j = self.original_coefficient_position(i)
            self.coefficients.append(original_coefficients[j].cpp_object())

        # Type checking coefficients
        if not all(isinstance(c, (cpp.function.GenericFunction))
                   for c in self.coefficients):
            coefficient_error = "Error while extracting coefficients. "
            raise TypeError(coefficient_error +
                            "Either provide a dict of cpp.function.GenericFunctions, " +
                            "or use Function to define your form.")

        for i in range(self.num_coefficients()):
            if isinstance(self.coefficients[i], cpp.function.GenericFunction):
                self.set_coefficient(i, self.coefficients[i])

        # Attach mesh (because function spaces and coefficients may be
        # empty lists)
        if not function_spaces:
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
