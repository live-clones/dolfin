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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <dolfin/common/Array.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/la/GenericVector.h>

namespace py = pybind11;

namespace expression_wrappers
{
  py::array_t<double> expression_eval(dolfin::Expression &instance,
                                      py::array_t<double> x)
  {
    // Create coordinate array
    const py::buffer_info x_info = x.request();
    const dolfin::Array<double> _x(x_info.shape[0], x.mutable_data());

    // Create object to hold data to be computed and returned
    py::array_t<double, py::array::c_style> values(instance.value_size());
    dolfin::Array<double> _values(instance.value_size(), values.mutable_data());

    // Evaluate basis
    instance.eval(_values, _x);
    return values;
  }

}


namespace dolfin_wrappers
{

  void function(py::module& m)
  {
    // Wrap dolfin::Expression
    py::class_<dolfin::Expression, std::shared_ptr<dolfin::Expression>>(m, "Expression")
      .def(py::init<std::size_t>())
      .def(py::init<std::size_t, std::size_t>())
      .def("eval", &expression_wrappers::expression_eval, "Evaluate Expression")
      .def("eval", (void (dolfin::Expression::*)(dolfin::Array<double>&, const dolfin::Array<double>&, const ufc::cell&) const) &dolfin::Expression::eval,
           "Evaluate Expression (cell version)")
      .def("test", []() { return "Expression test function"; });

    //-----------------------------------------------------------------------------
    // dolfin::Constant
    py::class_<dolfin::Constant, std::shared_ptr<dolfin::Constant>>(m, "Constant")
      .def(py::init<double>())
      .def(py::init<std::vector<double>>());

    //-----------------------------------------------------------------------------
    // dolfin::Function
    py::class_<dolfin::Function, std::shared_ptr<dolfin::Function>>(m, "Function")
      .def(py::init<std::shared_ptr<dolfin::FunctionSpace>>())
      .def("vector", (std::shared_ptr<dolfin::GenericVector> (dolfin::Function::*)())
           &dolfin::Function::vector);

    //-----------------------------------------------------------------------------
    // dolfin::FunctionSpace
    py::class_<dolfin::FunctionSpace, std::shared_ptr<dolfin::FunctionSpace>>(m, "FunctionSpace")
      .def(py::init<std::shared_ptr<dolfin::Mesh>, std::shared_ptr<dolfin::FiniteElement>,
           std::shared_ptr<dolfin::GenericDofMap>>());

  }

}
