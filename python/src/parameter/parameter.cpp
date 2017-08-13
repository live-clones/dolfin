// Copyright (C) 2017 Chris Richardson and Garth N. Wells
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
#include <pybind11/stl.h>

#include <dolfin/parameter/GlobalParameters.h>
#include <dolfin/parameter/Parameter.h>
#include <dolfin/parameter/Parameters.h>
#include "../mpi_interface.h"

namespace py = pybind11;

namespace dolfin_wrappers
{

  void parameter(py::module& m)
  {

    // dolfin::Parameters
    py::class_<dolfin::Parameters, std::shared_ptr<dolfin::Parameters>>
      (m, "Parameters")
      .def(py::init<>())
      .def(py::init<std::string>())
      .def(py::init<dolfin::Parameters>())
      // Use boost::variant to simplify
      .def("add", (void (dolfin::Parameters::*)(std::string, std::string)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(std::string, bool)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(std::string, int)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(std::string, double)) &dolfin::Parameters::add)
      .def("add", (void (dolfin::Parameters::*)(const dolfin::Parameters&)) &dolfin::Parameters::add)
      .def("get_range", [](dolfin::Parameters& self, std::string key)
           {
             const auto& p = self.find_parameter(key);
             std::set<std::string> range;
             p->get_range(range);
             return range;
           })
      .def("name", &dolfin::Parameters::name)
      .def("rename", &dolfin::Parameters::rename)
      .def("str", &dolfin::Parameters::str)
      .def("has_parameter", &dolfin::Parameters::has_parameter)
      .def("has_parameter_set", &dolfin::Parameters::has_parameter_set)
      .def("_update", &dolfin::Parameters::update)
      .def("_get_parameter", (dolfin::Parameter& (dolfin::Parameters::*)(std::string))
           &dolfin::Parameters::operator[], py::return_value_policy::reference)
      .def("_get_parameter_set", (dolfin::Parameters& (dolfin::Parameters::*)(std::string))
           &dolfin::Parameters::operator(), py::return_value_policy::reference)
      /*
        // With boost::variant, need to figure out how to set the
        // return policy per type in variant
      .def("__getitem__", [](dolfin::Parameters& self, std::string key)
           {
             mapbox::util::variant<bool, int, double, std::string, dolfin::Parameters*> v;
             auto param = self.find_parameter(key);
             if (auto param = self.find_parameter(key))
             {
               // FIXME: This will be simple once boost::variant is working
               //return param->value();
               auto _v = param->value();
               if (_v.which() == 1)
                 v = boost::get<bool>(_v);
               else if (_v.which() == 2)
                 v = boost::get<int>(_v);
               else if (_v.which() == 3)
                 v = boost::get<double>(_v);
               else if (_v.which() == 4)
                 v = boost::get<std::string>(_v);
               else
                 throw std::runtime_error("Cannot get parameter value");
               return v;
             }
             else if (auto param = self.find_parameter_set(key))
             {
               //v = &(*param);
               v = &self("key");
               return v;
             }
             else
               throw std::runtime_error("Invalid parameter key: " + key);
           }, py::return_value_policy::automatic_reference)
      */
      // FIXME: Implement checks and error handling below
      // FIXME: Can these functions be consolidated. Maybe boost::variant can help?
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, py::none value)
           {
             auto param = self.find_parameter(key);
             if (!param)
               throw std::runtime_error("Parameter not found in Parameters object");
             param->reset();
           }, "Reset Parameter (mark as unset) by setting to None.")
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, bool value)
           {
             auto param = self.find_parameter(key);
             *param = value;
           })
           }, py::arg(), py::arg().noconvert())
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, std::string value)
           {
             auto param = self.find_parameter(key);
             if (!param)
               throw std::runtime_error("Parameter not found in Parameters object");
             *param = value;
           })
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, int value)
           {
             auto param = self.find_parameter(key);
             *param = value;
           })
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, double value)
           {
             auto param = self.find_parameter(key);
             *param = value;
           })
      .def("__setitem__", [](dolfin::Parameters& self, std::string key, const dolfin::Parameters& other)
           {
             auto param = self.find_parameter_set(key);
             *param = other;
           })
      .def("parse", [](dolfin::Parameters& self, py::list argv)
           {
             if(argv.size() == 0)
               argv = py::module::import("sys").attr("argv").cast<py::list>();
             int argc = argv.size();
             std::vector<const char*> aptr;
             std::vector<std::string> a;
             for (auto q : argv)
             {
               a.push_back(q.cast<std::string>());
               aptr.push_back(a.back().c_str());
             }
             self.parse(argc, const_cast<char**>(aptr.data()));
           }, py::arg("argv")=py::list())
      .def("copy", [](dolfin::Parameters& self) { return dolfin::Parameters(self); })
      .def("assign", [](dolfin::Parameters& self, dolfin::Parameters& other) { self = other;});

    // dolfin::Parameter
    py::class_<dolfin::Parameter, std::shared_ptr<dolfin::Parameter>>(m, "Parameter")
      .def("value", &dolfin::Parameter::value)
      .def("__str__", &dolfin::Parameter::value_str);

    py::class_<dolfin::GlobalParameters, std::shared_ptr<dolfin::GlobalParameters>,
      dolfin::Parameters> (m, "GlobalParameters");

    // The global parameters (return a reference because there should
    // be only one instance)
    m.attr("parameters") = &dolfin::parameters;

  }

}
