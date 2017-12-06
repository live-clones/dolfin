"""Unit tests for the function library"""

# Copyright (C) 2017 Nate Sime
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

import pytest
from dolfin import *
from math import sin, cos, exp, tan
from numpy import array, zeros, float_
import numpy as np

from dolfin_utils.test import fixture, skip_in_parallel, skip_if_pybind11, skip_if_not_pybind11, has_pybind11


meshes = [UnitIntervalMesh(4), 
          UnitSquareMesh(4, 4, 'left/right'),
          UnitSquareMesh(4, 4, 'right/left'),
          UnitSquareMesh(4, 4, 'crossed'),
          UnitSquareMesh.create(4, 4, CellType.Type_quadrilateral),
          UnitCubeMesh(4, 4, 4),
          UnitCubeMesh.create(4, 4, 4, CellType.Type_hexahedron)]

if has_pybind11():
    from dolfin.mesh.meshfunction import _meshfunction_types
    meshfunction_types = set(_meshfunction_types.keys())
else:
    meshfunction_types = {"double", "size_t", "int", "bool"}

supported_dtypes = {"double"}
unsupported_dtypes = sorted(list(meshfunction_types - supported_dtypes))


@skip_if_not_pybind11
@pytest.mark.parametrize("mesh", meshes)
def test_mesh_expressions_cell(mesh):
    cf = MeshFunction("double", mesh, mesh.topology().dim(), 0.0)
    mvc = MeshValueCollection("double", mesh, mesh.topology().dim())

    for c in cells(mesh):
        cf[c] = c.midpoint()[0]
        mvc.set_value(c.index(), c.midpoint()[0])

    me_cf = MeshExpression(cf, degree=0)
    me_mvc = MeshExpression(mvc, degree=0)
    f = Expression("x[0]", degree=0)

    V = FunctionSpace(mesh, "CG", 1)
    u, v = TrialFunction(V), TestFunction(V)
    a = dot(grad(u), grad(v))*dx

    L_sym = f*v*dx
    L_me_cf = me_cf*v*dx
    L_me_mvc = me_mvc*v*dx

    bc = DirichletBC(V, Constant(0.0), "on_boundary")

    u_sym = Function(V)
    solve(a == L_sym, u_sym, bc)

    u_me_cf = Function(V)
    solve(a == L_me_cf, u_me_cf, bc)

    u_me_mvc = Function(V)
    solve(a == L_me_mvc, u_me_mvc, bc)

    for j in range(u_sym.vector().local_size()):
        assert(near(u_sym.vector()[j], u_me_cf.vector()[j], 1e-12))
        assert(near(u_sym.vector()[j], u_me_mvc.vector()[j], 1e-12))


@skip_if_not_pybind11
@pytest.mark.parametrize("mesh", meshes)
def test_mesh_expressions_facet(mesh):
    mesh.init_global(mesh.topology().dim() - 1)

    V = FunctionSpace(mesh, "CG", 1)
    u, v = TrialFunction(V), TestFunction(V)
    a = dot(grad(u), grad(v))*dx

    ff_d = MeshFunction("double", mesh, mesh.topology().dim()-1, 0.0)
    mvc = MeshValueCollection("double", mesh, mesh.topology().dim() - 1)
    ff_i = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

    ds = Measure("ds", subdomain_data=ff_i)
    L_sym = Constant(0.0)*v*ds

    exterior_facets_local = list(f.global_index() for f in facets(mesh) if f.exterior())
    proc_facets = mesh.mpi_comm().allgather(exterior_facets_local)
    exterior_facets = [facet for facets in proc_facets for facet in facets]
    num_exterior_facets = len(exterior_facets)

    for ext_facet_id in exterior_facets:
        L_sym += float(ext_facet_id + 1)*v*ds(ext_facet_id)

    for fa in facets(mesh):
        if fa.exterior():
            f_idx, f_val = fa.global_index(), float(fa.global_index()) + 1.0
            ff_d[fa] = f_val
            ff_i[fa] = f_idx

            c = Cell(mesh, fa.entities(mesh.topology().dim())[0])
            mvc.set_value(c.index(), c.index(fa), f_val)

    me_ff = MeshExpression(ff_d, degree=0)
    me_mvc = MeshExpression(mvc, degree=0)

    L_me_ff = me_ff*v*ds
    L_me_mvc = me_mvc*v*ds

    bc = DirichletBC(V, Constant(0.0), "near(x[0], 0.0)")

    u_sym = Function(V)
    solve(a == L_sym, u_sym, bc)

    u_me_ff = Function(V)
    solve(a == L_me_ff, u_me_ff, bc)

    u_me_mvc = Function(V)
    solve(a == L_me_mvc, u_me_mvc, bc)

    for j in range(u_sym.vector().local_size()):
        assert(near(u_sym.vector()[j], u_me_ff.vector()[j], 1e-12))
        assert(near(u_sym.vector()[j], u_me_mvc.vector()[j], 1e-12))


@skip_if_not_pybind11
@pytest.mark.parametrize("mesh", meshes)
@pytest.mark.parametrize("dtype", unsupported_dtypes)
def test_mesh_expressions_unsupported_type(mesh, dtype):

    ff = MeshFunction(dtype, mesh, mesh.topology().dim()-1, 0)
    cf = MeshFunction(dtype, mesh, mesh.topology().dim(), 0)

    with pytest.raises(TypeError):
        MeshExpression(cf)


@skip_if_not_pybind11
@pytest.mark.parametrize("mesh", meshes)
@pytest.mark.parametrize("dtype", supported_dtypes)
@pytest.mark.xfail
def test_mesh_expressions_unsupported_topology(mesh, dtype):
    if mesh.topology().dim() == 1:
        return

    for t_dim in range(mesh.topology().dim() - 1):
        mf = MeshFunction(dtype, mesh, t_dim)
        MeshExpression(mf)


@skip_if_not_pybind11
@pytest.mark.parametrize("mesh", meshes)
@pytest.mark.parametrize("dtype", supported_dtypes)
@pytest.mark.xfail
def test_mesh_expressions_mvc_unsupported_topology(mesh, dtype):
    if mesh.topology().dim() == 1:
        return

    for t_dim in range(mesh.topology().dim() - 1):
        mvc = MeshValueCollection(dtype, mesh, t_dim)
        MeshExpression(mvc)


@skip_if_not_pybind11
@pytest.mark.parametrize("mesh", meshes)
@pytest.mark.parametrize("dtype", supported_dtypes)
def test_mesh_expressions_unsupported_elements(mesh, dtype):
    if mesh.topology().dim() == 1:
        return

    for t_dim in range(mesh.topology().dim() - 1, mesh.topology().dim() + 1):
        mf = MeshFunction(dtype, mesh, t_dim)
        mvc = MeshValueCollection(dtype, mf)

        for data in [mf, mvc]:
            with pytest.raises(RuntimeError):
                MeshExpression(data, degree=1)

            with pytest.raises(RuntimeError):
                MeshExpression(data, element=FiniteElement("CG", mesh.ufl_cell(), 1))

            with pytest.raises(RuntimeError):
                MeshExpression(data, degree=1, element=FiniteElement("DG", mesh.ufl_cell(), 0))

            with pytest.raises(RuntimeError):
                MeshExpression(data, degree=0, element=FiniteElement("DG", mesh.ufl_cell(), 1))