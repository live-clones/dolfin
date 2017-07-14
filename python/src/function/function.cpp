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
#include <pybind11/eigen.h>


#include <dolfin/common/Array.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunction.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/la/GenericVector.h>

namespace py = pybind11;

namespace expression_wrappers
{
  /*
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
  */

}


namespace dolfin_wrappers
{

  void function(py::module& m)
  {
    // dolfin::GenericFunction
    py::class_<dolfin::GenericFunction, std::shared_ptr<dolfin::GenericFunction>>
      (m, "GenericFunction");

    // dolfin::MultiMeshFunction
    py::class_<dolfin::MultiMeshFunction, std::shared_ptr<dolfin::MultiMeshFunction>>
      (m, "MultiMeshFunction");

    //-------------------------------------------------------------------------

    m.def("make_dolfin_expression",
          [](std::uintptr_t e)
          {
            dolfin::Expression *p = reinterpret_cast<dolfin::Expression *>(e);
            return std::shared_ptr<const dolfin::Expression>(p);
          });

    // dolfin::Expression
    class PyExpression : public dolfin::Expression
    {
      using dolfin::Expression::Expression;

      void eval(Eigen::Ref<Eigen::VectorXd> values,
                const Eigen::Ref<Eigen::VectorXd> x) const override
      { PYBIND11_OVERLOAD(void, dolfin::Expression, eval, values, x); }

      void eval(Eigen::Ref<Eigen::VectorXd> values,
                const Eigen::Ref<Eigen::VectorXd> x,
                const ufc::cell& cell) const override
      { PYBIND11_OVERLOAD_NAME(void, dolfin::Expression, "eval_cell", eval, values, x, cell); }

    };


    py::class_<dolfin::Expression, PyExpression,
               std::shared_ptr<dolfin::Expression>, dolfin::GenericFunction,
               dolfin::Variable>
      (m, "Expression")
      .def(py::init<>())
      .def(py::init<std::size_t>())
      .def(py::init<std::size_t, std::size_t>())
      .def(py::init<std::vector<std::size_t>>())
      .def("eval", (void (dolfin::Expression::*)(Eigen::Ref<Eigen::VectorXd>, const Eigen::Ref<Eigen::VectorXd>) const)
           &dolfin::Expression::eval, "Evaluate Expression")
      .def("eval_cell", (void (dolfin::Expression::*)(Eigen::Ref<Eigen::VectorXd>, const Eigen::Ref<Eigen::VectorXd>, const ufc::cell&) const)
           &dolfin::Expression::eval,
           "Evaluate Expression (cell version)")
      .def("value_rank", &dolfin::Expression::value_rank)
      .def("value_dimension", &dolfin::Expression::value_dimension);


    //-----------------------------------------------------------------------------
    // dolfin::Constant
    py::class_<dolfin::Constant, std::shared_ptr<dolfin::Constant>, dolfin::Expression,
               dolfin::GenericFunction, dolfin::Variable>
      (m, "Constant")
      .def(py::init<double>())
      .def(py::init<std::vector<double>>());

    //-----------------------------------------------------------------------------
    // dolfin::Function
    py::class_<dolfin::Function, std::shared_ptr<dolfin::Function>, dolfin::GenericFunction>
      (m, "Function")
      .def(py::init<std::shared_ptr<dolfin::FunctionSpace>>())
      .def("__call__", [](dolfin::Function& self, std::vector<double>& p)
          {
            // FIXME - remove Array and replace with Eigen in DOLFIN
            const dolfin::Array<double> x(p.size(), p.data());
            Eigen::VectorXd values(self.value_size());
            dolfin::Array<double> _values(self.value_size(), values.data());
            self.eval(_values, x);
            return values;
          })
      .def("interpolate", &dolfin::Function::interpolate)
      .def("vector", (std::shared_ptr<dolfin::GenericVector> (dolfin::Function::*)())
           &dolfin::Function::vector);

    m.def("interpolate", [](const dolfin::GenericFunction& f,
                          std::shared_ptr<const dolfin::FunctionSpace> V)
          {
            auto g = std::make_shared<dolfin::Function>(V);
            g->interpolate(f);
            return g;
          });

    //-----------------------------------------------------------------------------
    // dolfin::FunctionSpace
    py::class_<dolfin::FunctionSpace, std::shared_ptr<dolfin::FunctionSpace>>
      (m, "FunctionSpace")
      .def(py::init<std::shared_ptr<dolfin::Mesh>, std::shared_ptr<dolfin::FiniteElement>,
           std::shared_ptr<dolfin::GenericDofMap>>());

  }

}
