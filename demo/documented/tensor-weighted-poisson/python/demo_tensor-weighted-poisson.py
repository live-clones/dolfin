"""This demo program solves Poisson's equation

    - div C grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 10*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0  for x = 0 or x = 1
du/dn(x, y) = 0  for y = 0 or y = 1

The conductivity C is a symmetric 2 x 2 matrix which
varies throughout the domain. In the left part of the
domain, the conductivity is

    C = ((1, 0.3), (0.3, 2))

and in the right part it is

    C = ((3, 0.5), (0.5, 4))

The data files where these values are stored are generated
by the program generate_data.py

This demo is dedicated to BF and Marius... ;-)
"""

# Copyright (C) 2009-2011 Anders Logg
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
# First added:  2009-12-16
# Last changed: 2011-06-28
# Begin demo

from dolfin import *

# Read mesh from file and create function space
mesh = Mesh("../unitsquare_32_32.xml.gz")
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Code for C++ evaluation of conductivity
conductivity_code = """

class Conductivity : public Expression
{
public:

  // Create expression with 3 components
  Conductivity() : Expression(3) {}

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& cell) const override
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.index;
    values[0] = (*c00)[cell_index];
    values[1] = (*c01)[cell_index];
    values[2] = (*c11)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<MeshFunction<double>> c00;
  std::shared_ptr<MeshFunction<double>> c01;
  std::shared_ptr<MeshFunction<double>> c11;

};
"""

conductivity_pybind11 = """
  py::class_<dolfin::Conductivity, std::shared_ptr<dolfin::Conductivity>, dolfin::Expression>
    (m, "Conductivity", py::dynamic_attr())
    .def(py::init<>())
    .def("eval", &dolfin::Conductivity::eval)
    .def_readwrite("c00", &dolfin::Conductivity::c00)
    .def_readwrite("c01", &dolfin::Conductivity::c01)
    .def_readwrite("c11", &dolfin::Conductivity::c11);
"""

conductivity_class_data = {'cpp_code': conductivity_code, 'pybind11_code': conductivity_pybind11}

# Define conductivity expression and matrix
c00 = MeshFunction("double", mesh, "../unitsquare_32_32_c00.xml.gz")
c01 = MeshFunction("double", mesh, "../unitsquare_32_32_c01.xml.gz")
c11 = MeshFunction("double", mesh, "../unitsquare_32_32_c11.xml.gz")

class UserConductivity(UserExpression):
    def value_shape(self):
        return (3,)

c = UserConductivity(degree=0)
cc = compile_cpp_code(conductivity_class_data).Conductivity()
cc.c00 = c00
cc.c01 = c01
cc.c11 = c11
c._cpp_object = cc

C = as_matrix(((c[0], c[1]), (c[1], c[2])))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
a = inner(C*grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in VTK format
file = File("poisson.pvd")
file << u

# Plot solution
import matplotlib.pyplot as plt
plot(u)
plt.show()
