// Copyright (C) 2017 Chris Richardson
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

#include <dolfin/parameter/GlobalParameters.h>
#include <dolfin/parameter/Parameter.h>
#include <dolfin/parameter/Parameters.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void parameter(py::module& m)
  {

    py::class_<dolfin::Parameters, std::shared_ptr<dolfin::Parameters>>
      (m, "Parameters")
      .def(py::init<>())
      .def(py::init<std::string>())
      .def(py::init<dolfin::Parameters>())
      .def("add", (void (dolfin::Parameters::*)(std::string, std::string)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(std::string, int)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(std::string, double)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(const dolfin::Parameters&)) &dolfin::Parameters::add)
      .def("name", &dolfin::Parameters::name)
      .def("rename", &dolfin::Parameters::rename)
      .def("update", &dolfin::Parameters::update)
      .def("has_parameter", &dolfin::Parameters::has_parameter)
      .def("has_parameter_set", &dolfin::Parameters::has_parameter_set)
      .def("_get_parameter", (dolfin::Parameter& (dolfin::Parameters::*)(std::string))
           &dolfin::Parameters::operator[], py::return_value_policy::reference)
      .def("_get_parameter_set", (dolfin::Parameters& (dolfin::Parameters::*)(std::string))
           &dolfin::Parameters::operator(), py::return_value_policy::reference)
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, std::string value)
           {
             auto param = self.find_parameter(key);
             *param = value;
           })
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, bool value)
           {
             auto param = self.find_parameter(key);
             *param = value;
           })
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, int value)
           {
             auto param = self.find_parameter(key);
             *param = value;
           })
      .def("copy", [](dolfin::Parameters& self) { return dolfin::Parameters(self); })
      .def("assign", [](dolfin::Parameters& self, dolfin::Parameters& other) { self = other;}) ;


    py::class_<dolfin::Parameter, std::shared_ptr<dolfin::Parameter>>
      (m, "Parameter");

    py::class_<dolfin::IntParameter, std::shared_ptr<dolfin::IntParameter>,
      dolfin::Parameter>
      (m, "IntParameter")
      .def("value", [](dolfin::IntParameter& self) { return int(self); })
      .def("__str__", &dolfin::IntParameter::value_str);

    py::class_<dolfin::DoubleParameter, std::shared_ptr<dolfin::DoubleParameter>,
      dolfin::Parameter>
      (m, "DoubleParameter")
      .def("value", [](dolfin::DoubleParameter& self) { return double(self); })
      .def("__str__", &dolfin::DoubleParameter::value_str);

    py::class_<dolfin::StringParameter, std::shared_ptr<dolfin::StringParameter>,
               dolfin::Parameter>
      (m, "StringParameter")
      .def("value", [](dolfin::StringParameter& self) { return std::string(self); })
      .def("__str__", &dolfin::StringParameter::value_str);

    py::class_<dolfin::BoolParameter, std::shared_ptr<dolfin::BoolParameter>,
      dolfin::Parameter>
      (m, "BoolParameter")
      .def("value", [](dolfin::BoolParameter& self) { return bool(self); })
      .def("__str__", &dolfin::BoolParameter::value_str);

    py::class_<dolfin::GlobalParameters, std::shared_ptr<dolfin::GlobalParameters>,
      dolfin::Parameters> (m, "GlobalParameters");

    m.attr("parameters") = dolfin::parameters;

  }

}
