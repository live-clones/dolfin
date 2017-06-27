// Copyright (C) 2017 Garth N. Wells
//
// This file is part of FFC.
//
// FFC is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FFC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FFC. If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <dolfin/mesh/Mesh.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void mesh(py::module& m)
  {
    // Wrap dolfin::Mesh class
    py::class_<dolfin::Mesh, std::shared_ptr<dolfin::Mesh>> Mesh(m, "Mesh", "DOLFIN Mesh object");

    // Constructors
    Mesh.def(py::init<>());

    // Mesh member functions
    Mesh.def("num_entities", &dolfin::Mesh::num_entities, "Number of mesh entities");
  }

}
