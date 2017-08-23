# -*- coding: utf-8 -*-

import dolfin.cpp as cpp
from . import get_pybind_include
from dolfin.function.expression import BaseExpression, _select_element

from six import string_types
import hashlib
import dijitso

def jit_generate(class_data, module_name, signature, parameters):

#FIXME: decide what needs to be in boilerplate code, if anything.

    template_code = """

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin.h>
#include <Eigen/Dense>

namespace dolfin
{{
   {cpp_code}
}}

PYBIND11_MODULE({signature}, m)
{{
   {pybind11_code}
}}
"""

    code_c = template_code.format(cpp_code=class_data["cpp_code"],
                                  pybind11_code=class_data["pybind11_code"],
                                  signature=signature)
    code_h = ""
    depends = []

    return code_h, code_c, depends

def compile_cpp_code(class_data):
    """Compile a user C(++) string to a Python object with pybind11.
       Note this is still experimental."""

    import pkgconfig
    if not pkgconfig.exists('dolfin'):
        raise RuntimeError("Could not find DOLFIN pkg-config file. Please make sure appropriate paths are set.")

    # Get pkg-config data
    d = pkgconfig.parse('dolfin')

    # Set compiler/build options
    # FIXME: need to locate Python libs and pybind11 - this works on my Mac...
    from distutils import sysconfig
    params = dijitso.params.default_params()
    pyversion = "python" + sysconfig.get_config_var("LDVERSION")
    params['cache']['lib_prefix'] = ""
    params['cache']['lib_basename'] = ""
    params['cache']['lib_loader'] = "import"
    params['build']['include_dirs'] = d["include_dirs"] + get_pybind_include() + [sysconfig.get_config_var("INCLUDEDIR") + "/" + pyversion]
    params['build']['libs'] = d["libraries"] + [ pyversion ]
    params['build']['lib_dirs'] = d["library_dirs"] + [sysconfig.get_config_var("LIBDIR")]
    params['build']['cxxflags'] += ('-fno-lto',)

    # FIXME: should make this generalised to any library
    if cpp.common.has_petsc():
        import os
        params['build']['cxxflags'] += ('-DHAS_PETSC',)
        params['build']['libs'] += ['petsc']
        params['build']['lib_dirs'] += [os.environ["PETSC_DIR"] + "/lib"]

    print(params)

    module_hash = hashlib.md5((class_data["cpp_code"] + class_data["pybind11_code"]).encode('utf-8')).hexdigest()
    module_name = "dolfin_cpp_module_" + module_hash

    module, signature = dijitso.jit(class_data,
                                    module_name, params,
                                    generate=jit_generate)

    print('pybind11 JIT: ', module)

    print(dir(module))

    return module

class CompiledExpressionPyBind11(BaseExpression):
    """JIT Expressions"""

    def __init__(self, cpp_code=None, *args, **kwargs):

        # Remove arguments that are used in Expression creation
        element = kwargs.pop("element", None)
        degree = kwargs.pop("degree", None)
        cell = kwargs.pop("cell", None)
        domain = kwargs.pop("domain", None)
        name = kwargs.pop("name", None)
        label = kwargs.pop("label", None)
        mpi_comm = kwargs.pop("mpi_comm", None)

        # Save properties for checking later
        self._properties = kwargs
        for k in self._properties:
            if not isinstance(k, string_types):
                raise KeyError("Invalid key")

        if cpp_code is not None:
            self._cpp_object = compile_expression(cpp_code, self._properties)

        # Deduce element type if not provided
        if element is None:
            if degree is None:
                raise KeyError("Must supply element or degree")
            value_shape = tuple(self.value_dimension(i)
                                for i in range(self.value_rank()))
            element = _select_element(family=None, cell=None, degree=degree,
                                      value_shape=value_shape)

        BaseExpression.__init__(self, cell=cell, element=element, domain=domain,
                                name=name, label=label)

    # This is added dynamically in the intialiser to allow checking of
    # eval in user classes.
    def __getattr__(self, name):
        "Pass attributes through to (JIT compiled) Expression object"
        if name in self._properties.keys():
            return getattr(self._cpp_object, name)
        else:
            raise(AttributeError)

    def __setattr__(self, name, value):
        if name.startswith("_"):
            super().__setattr__(name, value)
        elif name in self._properties.keys():
            setattr(self._cpp_object, name, value)
