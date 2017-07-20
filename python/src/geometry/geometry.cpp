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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>

#include <dolfin/geometry/intersect.h>
#include <dolfin/geometry/BoundingBoxTree.h>
#include <dolfin/geometry/MeshPointIntersection.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/mesh/Mesh.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void geometry(py::module& m)
  {
    // dolfin::BoundingBoxTree
    py::class_<dolfin::BoundingBoxTree, std::shared_ptr<dolfin::BoundingBoxTree>>
      (m, "BoundingBoxTree")
      .def(py::init<>())
      .def("build", (void (dolfin::BoundingBoxTree::*)(const dolfin::Mesh&))
           &dolfin::BoundingBoxTree::build)
      .def("build", (void (dolfin::BoundingBoxTree::*)(const dolfin::Mesh&, std::size_t))
           &dolfin::BoundingBoxTree::build)
      .def("compute_collisions", (std::vector<unsigned int> (dolfin::BoundingBoxTree::*)(const dolfin::Point&) const)
           &dolfin::BoundingBoxTree::compute_collisions)
      .def("compute_collisions",
           (std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
            (dolfin::BoundingBoxTree::*)(const dolfin::BoundingBoxTree&) const)
           &dolfin::BoundingBoxTree::compute_collisions)
      .def("compute_entity_collisions", (std::vector<unsigned int> (dolfin::BoundingBoxTree::*)(const dolfin::Point&) const)
           &dolfin::BoundingBoxTree::compute_entity_collisions)
      .def("compute_entity_collisions",
           (std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
            (dolfin::BoundingBoxTree::*)(const dolfin::BoundingBoxTree&) const)
            &dolfin::BoundingBoxTree::compute_entity_collisions)
      .def("compute_first_collision", &dolfin::BoundingBoxTree::compute_first_collision)
      .def("compute_first_entity_collision", &dolfin::BoundingBoxTree::compute_first_entity_collision)
      .def("compute_closest_entity", &dolfin::BoundingBoxTree::compute_closest_entity);


    // dolfin::Point
    py::class_<dolfin::Point>(m, "Point")
      .def(py::init<>())
      .def(py::init<double, double, double>())
      .def(py::init<double, double>())
      .def(py::init<double>())
      .def("__init__",
           [](dolfin::Point& instance, py::array_t<double> x)
           {
             auto b = x.request();
             assert(b.shape.size() == 1);
             assert(b.shape.size()[0] <= 3);
             new (&instance) dolfin::Point(b.shape[0], x.data());
           })
      .def("__getitem__", [](const dolfin::Point& self, std::size_t index)
           { if (index > 2)
               throw py::index_error("Out of range");
             return self[index]; })
      .def("__setitem__", [](dolfin::Point& self, std::size_t index, double value)
           { if (index > 2)
               throw py::index_error("Out of range");
             self[index] = value; })
      .def("__add__", [](const dolfin::Point& self, const dolfin::Point& other)
           { return self + other; })
      .def("__sub__", [](const dolfin::Point& self, const dolfin::Point& other)
           { return self - other; })
      .def("array", [](dolfin::Point& self)
           { return Eigen::Map<Eigen::Vector3d>(self.coordinates()); })
      .def("norm", &dolfin::Point::norm)
      .def("x", &dolfin::Point::x)
      .def("y", &dolfin::Point::y)
      .def("z", &dolfin::Point::z)
      .def("distance", &dolfin::Point::distance);

    // dolfin::MeshPointIntersection
    py::class_<dolfin::MeshPointIntersection, std::shared_ptr<dolfin::MeshPointIntersection>>
      (m, "MeshPointIntersection")
      .def("intersected_cells", &dolfin::MeshPointIntersection::intersected_cells);

    m.def("intersect", &dolfin::intersect);

  }
}
