// Copyright (C) 2017 Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#include <dolfin/geometry/Point.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void geometry(py::module& m)
  {
    // Wrap dolfin::Point
    py::class_<dolfin::Point>(m, "Point")
      .def(py::init<>())
      .def(py::init<double, double, double>())
      .def(py::init<double, double>())
      .def(py::init<double>())
      .def("__getitem__", [](const dolfin::Point& self, std::size_t index)
           { return self[index]; })
      .def("__setitem__", [](dolfin::Point& self, std::size_t index, double value)
           { self[index] = value; })
      .def("__add__", [](const dolfin::Point& self, const dolfin::Point& other)
           { return self+other; })
      .def("__sub__", [](const dolfin::Point& self, const dolfin::Point& other)
           { return self-other; })
      .def("array", [](dolfin::Point& self)
           { return Eigen::Map<Eigen::Vector3d>(self.coordinates()); })
      .def("norm", &dolfin::Point::norm);
  }

}
