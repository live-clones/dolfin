# -*- coding: utf-8 -*-
"""This module provides functionality for form assembly in Python,
corresponding to the C++ assembly and PDE classes.

The C++ :py:class:`assemble <dolfin.cpp.assemble>` function
(renamed to cpp_assemble) is wrapped with an additional
preprocessing step where code is generated using the
FFC JIT compiler.

The C++ PDE classes are reimplemented in Python since the C++ classes
rely on the dolfin::Form class which is not used on the Python side.

"""

# Copyright (C) 2007-2015 Anders Logg
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

import ufl_legacy as ufl
import dolfin.cpp as cpp
from dolfin.fem.form import Form
from dolfin import MPI
from dolfin.function.multimeshfunction import MultiMeshFunction

from ufl_legacy.form import sub_forms_by_domain

__all__ = ["assemble", "assemble_mixed", "assemble_local", "assemble_system",
           "assemble_multimesh", "SystemAssembler"]


def _create_dolfin_form(form, form_compiler_parameters=None,
                        function_spaces=None):
    if form is None:
        return form
    # First check if we got a cpp.Form
    elif isinstance(form, cpp.fem.Form):

        # Check that jit compilation has already happened
        if not hasattr(form, "_compiled_form"):
            raise TypeError("Expected a dolfin form to have a _compiled_form attribute.")

        # Warn that we don't use the parameters if we get any
        if form_compiler_parameters is not None:
            cpp.warning("Ignoring form_compiler_parameters when passed a dolfin Form!")
        return form
    elif isinstance(form, ufl.Form):
        return Form(form,
                    form_compiler_parameters=form_compiler_parameters,
                    function_spaces=function_spaces)
    else:
        raise TypeError("Invalid form type %s" % (type(form),))


def assemble_local(form, cell, form_compiler_parameters=None):
    """JIT assemble_local"""
    # Create dolfin Form object
    if isinstance(form, cpp.fem.Form):
        dolfin_form = form
    else:
        dolfin_form = _create_dolfin_form(form, form_compiler_parameters)
    result = cpp.fem.assemble_local(dolfin_form, cell)
    if result.shape[1] == 1:
        if result.shape[0] == 1:
            result = result[0][0]
        else:
            result = result.reshape((result.shape[0]))
    return result


def assemble(form, tensor=None, form_compiler_parameters=None,
             add_values=False, finalize_tensor=True,
             keep_diagonal=False, backend=None):
    """Assemble the given form and return the corresponding tensor.

    *Arguments*
        Depending on the input form, which may be a functional, linear
        form, bilinear form or higher rank form, a scalar value, a vector,
        a matrix or a higher rank tensor is returned.

    In the simplest case, no additional arguments are needed. However,
    additional arguments may and must in some cases be provided as
    outlined below.

    The ``form`` can be either a UFL Form or a DOLFIN Form.

    If the form defines integrals over different subdomains,
    :py:class:`MeshFunctions <dolfin.cpp.MeshFunction>` over the
    corresponding topological entities defining the subdomains can be
    provided.

    The implementation of the returned tensor is determined by the
    default linear algebra backend. This can be overridden by
    specifying a different backend.

    Each call to assemble() will create a new tensor. If the
    ``tensor`` argument is provided, this will be used instead.
    Sparsity pattern of ``tensor`` is reset iff ``tensor.empty()``
    holds.

    If the ``keep_diagonal`` is set to True, assembler ensures that
    potential zeros on a matrix diagonal are kept in sparsity pattern
    so every diagonal entry can be changed in a future (for example by
    ident() or ident_zeros()).

    Specific form compiler parameters can be provided by the
    ``form_compiler_parameters`` argument. Form compiler parameters
    can also be controlled using the global parameters stored in
    parameters["form_compiler"].

    *Examples of usage*
        The standard stiffness matrix ``A`` and a load vector ``b``
        can be assembled as follows:

        .. code-block:: python

            A = assemble(inner(grad(u),grad(v))*dx)
            b = assemble(f*v*dx)

        To prescribe the domains over which integrals will be
        evaluated, create a Measure with the MeshFunction passed as
        subdomain_data.  For instance, using a mesh function marking
        parts of the boundary:

        .. code-block:: python

            # MeshFunction marking boundary parts
            boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
            # ... fill values in boundary_markers

            # Measures with references to cell and boundary markers
            ds = Measure("ds", subdomain_data=boundary_markers)

            # Sample variational forms
            a = inner(grad(u), grad(v))*dx + p*u*v*ds(0)
            L = f*v*dx - g*v*ds(1) + p*q*v*ds(0)

            A = assemble(a)
            b = assemble(L)

        For interior facet integrals (dS), cell markers can be used to
        specify which cell is '+' and which is '-'. The '+' and '-'
        sides are chosen such that the cell marker value in the cell
        at the '+' side cell is larger than the cell marker value in
        the cell at the '-' side cell. If the values are equal or the
        cell markers are not provided, the sides are chosen
        arbitrarily.

        Note that currently, cell markers must be associated with a
        cell type integral (dx), and if you don't have such an
        integral a workaround is to add the integral of something over
        an empty domain such as 'f*dx(99)' with 99 a number not
        occuring among the cell markers. A better interface to handle
        this case will be provided later.

        .. code-block:: python

            # MeshFunctions marking boundary and cell parts
            boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
            cell_markers = MeshFunction("size_t", mesh, mesh.topology().dim())
            # ... fill values in boundary_markers

            # Measures with references to cell and boundary markers
            ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
            dx = Measure("dx", domain=mesh, subdomain_data=cell_markers)

            # Sample variational forms
            a = inner(grad(u), grad(v))*dx + p*u*v*ds(0)
            L = v*dx(99) - g*v*ds(1) + p*q*v*ds(0)

            A = assemble(a)
            b = assemble(L)

        To ensure that the assembled matrix has the right type, one may use
        the ``tensor`` argument:

        .. code-block:: python

            A = PETScMatrix()
            assemble(a, tensor=A)

        The form ``a`` is now assembled into the PETScMatrix ``A``.

    """

    # Create dolfin Form object
    if isinstance(form, cpp.fem.Form):
        dolfin_form = form
    else:
        dolfin_form = _create_dolfin_form(form, form_compiler_parameters)

    # Create tensor
    comm = dolfin_form.mesh().mpi_comm()
    tensor = _create_tensor(comm, form, dolfin_form.rank(), backend, tensor)

    # Create C++ assembler
    assembler = cpp.fem.Assembler()

    # Set assembler options
    assembler.add_values = add_values
    assembler.finalize_tensor = finalize_tensor
    assembler.keep_diagonal = keep_diagonal

    # Call C++ assemble
    assembler.assemble(tensor, dolfin_form)

    # Convert to float for scalars
    if dolfin_form.rank() == 0:
        tensor = tensor.get_scalar_value()

    # Return value
    return tensor

# JIT assembler


def assemble_mixed(form,
                   tensor=None,
                   form_compiler_parameters=None,
                   add_values=False,
                   finalize_tensor=True,
                   keep_diagonal=False,
                   backend=None):

    # Create dolfin Form object referencing all data needed by assembler
    if isinstance(form, cpp.fem.Form):
        dolfin_forms = [form]
    else:
        dolfin_forms = []
        for subform in sub_forms_by_domain(form):
            dolfin_forms.append(_create_dolfin_form(subform, form_compiler_parameters))

    # Create tensor
    comm = dolfin_forms[0].mesh().mpi_comm()
    tensor = _create_tensor(comm, form, dolfin_forms[0].rank(), backend, tensor)

    # Create C++ mixed assembler
    assembler = cpp.fem.MixedAssembler()

    # Set assembler options
    assembler.add_values = add_values
    assembler.finalize_tensor = finalize_tensor
    assembler.keep_diagonal = keep_diagonal

    # Call C++ assemble
    for k, dolfin_form in enumerate(dolfin_forms):
        assembler.add_values = bool(k > 0)
        assembler.assemble(tensor, dolfin_form)

    # Convert to float for scalars
    if dolfin_forms[0].rank() == 0:
        tensor = tensor.get_scalar_value()

    # Return value
    return tensor

# JIT multimesh assembler


def assemble_multimesh(form,
                       tensor=None,
                       form_compiler_parameters=None,
                       backend=None):
    "Assemble the given multimesh form and return the corresponding tensor."

    # The form that comes in is (by construction in function.Argument)
    # defined on the first part of the multimesh. We now need to create
    # the DOLFIN Forms with the proper function spaces for each part.

    # FIXME: This code makes a number of assumptions and will need to
    # be revisited and improved.

    # Warn that we don't use the parameters if we get any
    if form_compiler_parameters is not None:
        cpp.warning("Ignoring form_compiler_parameters when passed a dolfin Form!")
    form_compiler_parameters = None

    # Extract arguments and multimesh function space
    coefficients = form.coefficients()
    arguments = form.arguments()

    # Extract rank
    rank = len(arguments)

    # Extract multimesh function spaces for arguments
    V_multi = [v._V_multi for v in arguments]

    # Extract number of parts, the multimesh and create the multimesh form
    num_parts = None
    if rank > 0:
        num_parts = V_multi[0].num_parts()
        multimesh_form = cpp.fem.MultiMeshForm(*V_multi)
        multimesh = V_multi[0].multimesh()
    elif len(coefficients) > 0:
        for coeff in coefficients:
            # Only create these variables once
            if isinstance(coeff, MultiMeshFunction):
                multimesh = coeff.function_space().multimesh()
                num_parts = coeff.function_space().num_parts()
                multimesh_form = cpp.fem.MultiMeshForm(multimesh)
                break

    if not num_parts:
        # Handle the case Constant(1)*dx(domain=multimesh)
        multimesh = form.ufl_domains()[0].ufl_cargo()
        num_parts = multimesh.num_parts()
        multimesh_form = cpp.fem.MultiMeshForm(multimesh)

    # Build multimesh DOLFIN form
    for part in range(num_parts):
        # Extract standard function spaces for all arguments on
        # current part
        function_spaces = [V_multi[i].part(part) for i in range(rank)]
        # Wrap standard form
        dolfin_form = _create_dolfin_form(form,
                                          form_compiler_parameters,
                                          function_spaces)

        # Setting coefficients for the multimesh form
        for i in range(len(coefficients)):
            if isinstance(coefficients[i], MultiMeshFunction):
                coeff = coefficients[i].part(part)
            else:
                coeff = coefficients[i]
            # Developer note: This may be done more elegantly by modifiying
            # _create_dolfin_form
            dolfin_form.set_coefficient(i, coeff._cpp_object)
            dolfin_form.coefficients[i] = coeff

        # Add standard mesh to the standard form and the
        # standard form to the multimesh form
        dolfin_form.set_mesh(multimesh.part(part))
        multimesh_form.add(dolfin_form)

    for i, coeff in enumerate(coefficients):
        if isinstance(coeff, MultiMeshFunction):
            multimesh_form.set_multimesh_coefficient(i, coeff._cpp_object)

    # Build multimesh form
    multimesh_form.build()

    # Create tensor
    comm = MPI.comm_world
    tensor = _create_tensor(comm, form, rank, backend, tensor)

    # Call C++ assemble function
    assembler = cpp.fem.MultiMeshAssembler()
    assembler.assemble(tensor, multimesh_form)

    # Convert to float for scalars
    if rank == 0:
        tensor = tensor.get_scalar_value()

    # Return value
    return tensor


def assemble_system(A_form, b_form, bcs=None, x0=None,
                    form_compiler_parameters=None, add_values=False,
                    finalize_tensor=True, keep_diagonal=False,
                    A_tensor=None, b_tensor=None, backend=None):
    """Assemble form(s) and apply any given boundary conditions in a
    symmetric fashion and return tensor(s).

    The standard application of boundary conditions does not
    necessarily preserve the symmetry of the assembled matrix. In
    order to perserve symmetry in a system of equations with boundary
    conditions, one may use the alternative assemble_system instead of
    multiple calls to :py:func:`assemble
    <dolfin.fem.assembling.assemble>`.

    *Examples of usage*

       For instance, the statements

       .. code-block:: python

           A = assemble(a)
           b = assemble(L)
           bc.apply(A, b)

       can alternatively be carried out by

       .. code-block:: python

           A, b = assemble_system(a, L, bc)

       The statement above is valid even if ``bc`` is a list of
       :py:class:`DirichletBC <dolfin.fem.bcs.DirichletBC>`
       instances. For more info and options, see :py:func:`assemble
       <dolfin.fem.assembling.assemble>`.

    """
    # Create dolfin Form objects referencing all data needed by
    # assembler
    A_dolfin_form = _create_dolfin_form(A_form, form_compiler_parameters)
    b_dolfin_form = _create_dolfin_form(b_form, form_compiler_parameters)

    # Create tensors
    comm_A = A_dolfin_form.mesh().mpi_comm()
    comm_b = A_dolfin_form.mesh().mpi_comm()
    A_tensor = _create_tensor(comm_A, A_form, A_dolfin_form.rank(), backend,
                              A_tensor)
    b_tensor = _create_tensor(comm_b, b_form, b_dolfin_form.rank(), backend,
                              b_tensor)

    # Check bcs
    bcs = _wrap_in_list(bcs, 'bcs', cpp.fem.DirichletBC)

    # Call C++ assemble function
    assembler = cpp.fem.SystemAssembler(A_dolfin_form, b_dolfin_form, bcs)
    assembler.add_values = add_values
    assembler.finalize_tensor = finalize_tensor
    assembler.keep_diagonal = keep_diagonal
    if x0 is not None:
        assembler.assemble(A_tensor, b_tensor, x0)
    else:
        assembler.assemble(A_tensor, b_tensor)

    return A_tensor, b_tensor


def _wrap_in_list(obj, name, types=type):
    if obj is None:
        lst = []
    elif hasattr(obj, '__iter__'):
        lst = list(obj)
    else:
        lst = [obj]
    for obj in lst:
        if not isinstance(obj, types):
            raise TypeError("expected a (list of) %s as '%s' argument" %
                            (str(types), name))
    return lst


def _create_tensor(mpi_comm, form, rank, backend, tensor):
    """Create tensor for form"""

    # Check if tensor is supplied by user
    if tensor is not None:
        return tensor

    # Check backend argument
    if (backend is not None) and (not isinstance(backend, cpp.la.GenericLinearAlgebraFactory)):
        raise TypeError("Provide a GenericLinearAlgebraFactory as 'backend'")

    # Create tensor
    if rank == 0:
        tensor = cpp.la.Scalar(mpi_comm)
    elif rank == 1:
        if backend:
            tensor = backend.create_vector(mpi_comm)
        else:
            tensor = cpp.la.Vector(mpi_comm)
    elif rank == 2:
        if backend:
            tensor = backend.create_matrix(mpi_comm)
        else:
            tensor = cpp.la.Matrix(mpi_comm)
    else:
        raise RuntimeError("Unable to create tensors of rank %d." % rank)

    return tensor


class SystemAssembler(cpp.fem.SystemAssembler):
    __doc__ = cpp.fem.SystemAssembler.__doc__

    def __init__(self, A_form, b_form, bcs=None,
                 form_compiler_parameters=None):
        """
        Create a SystemAssembler

        * Arguments *
           a (ufl.Form, _Form_)
              Bilinear form
           L (ufl.Form, _Form_)
              Linear form
           bcs (_DirichletBC_)
              A list or a single DirichletBC (optional)
        """

        if isinstance(A_form, list) and isinstance(b_form, list):
            A_dolfin_forms = [_create_dolfin_form(f, form_compiler_parameters) for f in A_form]
            b_dolfin_forms = [_create_dolfin_form(f, form_compiler_parameters) for f in b_form]

            # Call C++ assemble function
            cpp.fem.SystemAssembler.__init__(self, A_dolfin_forms, b_dolfin_forms,
                                             bcs)

            # Keep Python counterpart of bcs (and Python object it owns)
            # alive
            self._bcs = bcs
        else:
            # Create dolfin Form objects referencing all data needed by
            # assembler
            A_dolfin_form = _create_dolfin_form(A_form, form_compiler_parameters)
            b_dolfin_form = _create_dolfin_form(b_form, form_compiler_parameters)

            # Check bcs
            bcs = _wrap_in_list(bcs, 'bcs', cpp.fem.DirichletBC)

            # Call C++ assemble function
            cpp.fem.SystemAssembler.__init__(self, A_dolfin_form, b_dolfin_form,
                                             bcs)

            # Keep Python counterpart of bcs (and Python object it owns)
            # alive
            self._bcs = bcs
