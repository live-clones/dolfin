/// Copyright (C) 2017 Garth N. Wells
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
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionAXPY.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunction.h>
#include <dolfin/function/SpecialFunctions.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Mesh.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void function(py::module& m)
  {
    // Delcare ufc::cell
    py::class_<ufc::cell, std::shared_ptr<ufc::cell>>(m, "ufc_cell")
      .def_readonly("index", &ufc::cell::index);

    // Delcare ufc::function
    py::class_<ufc::function, std::shared_ptr<ufc::function>>(m, "ufc_function");

    // dolfin::GenericFunction
    py::class_<dolfin::GenericFunction, std::shared_ptr<dolfin::GenericFunction>,
               ufc::function, dolfin::Variable>
      (m, "GenericFunction")
      .def("value_dimension", &dolfin::GenericFunction::value_dimension)
      .def("eval", (void (dolfin::GenericFunction::*)(Eigen::Ref<Eigen::VectorXd>,
                                                      const Eigen::Ref<Eigen::VectorXd>, const ufc::cell&) const)
           &dolfin::GenericFunction::eval,
           "Evaluate GenericFunction (cell version)")
      .def("eval", (void (dolfin::GenericFunction::*)(Eigen::Ref<Eigen::VectorXd>, const Eigen::Ref<Eigen::VectorXd>) const)
           &dolfin::GenericFunction::eval, py::arg("values"), py::arg("x"), "Evaluate GenericFunction")
      .def("compute_vertex_values", [](dolfin::GenericFunction& self, const dolfin::Mesh& mesh)
           { std::vector<double> values;
             self.compute_vertex_values(values, mesh);
             return py::array_t<double>(values.size(), values.data());
           })
      .def("function_space", &dolfin::GenericFunction::function_space);

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

    py::class_<dolfin::Expression, PyExpression, std::shared_ptr<dolfin::Expression>,
               dolfin::GenericFunction>(m, "Expression")
      .def(py::init<>())
      .def(py::init<std::size_t>())
      .def(py::init<std::size_t, std::size_t>())
      .def(py::init<std::vector<std::size_t>>())
      .def("__call__", [](const dolfin::Expression& self, const Eigen::Ref<Eigen::VectorXd> x)
           {
             Eigen::VectorXd f(self.value_size());
             self.eval(f, x);
             return f;
           })
      .def("__call__", [](const dolfin::Expression& self, const dolfin::Point& p)
           {
             const Eigen::Map<Eigen::VectorXd> x(const_cast<double*>(p.coordinates()), 3);
             Eigen::VectorXd f(self.value_size());
             self.eval(f, x);
             return f;
           })
      .def("eval", (void (dolfin::Expression::*)(Eigen::Ref<Eigen::VectorXd>,
                                                      const Eigen::Ref<Eigen::VectorXd>, const ufc::cell&) const)
           &dolfin::Expression::eval,
           "Evaluate Expression (cell version)")
      .def("eval", (void (dolfin::Expression::*)(Eigen::Ref<Eigen::VectorXd>, const Eigen::Ref<Eigen::VectorXd>) const)
           &dolfin::Expression::eval, py::arg("values"), py::arg("x"), "Evaluate Expression")
      .def("value_rank", &dolfin::Expression::value_rank)
      .def("value_dimension", &dolfin::Expression::value_dimension)
      .def("get_property", &dolfin::Expression::get_property)
      .def("set_property", &dolfin::Expression::set_property);


    //-----------------------------------------------------------------------------
    // dolfin::Constant
    py::class_<dolfin::Constant, std::shared_ptr<dolfin::Constant>, dolfin::Expression>
      (m, "Constant")
      .def(py::init<double>())
      .def(py::init<std::vector<double>>())
      .def(py::init<std::vector<std::size_t>, std::vector<double>>())
      .def("values", [](const dolfin::Constant& self)
           { auto v =  self.values(); return py::array_t<double>(v.size(), v.data()); })
      .def("__float__", [](const dolfin::Constant& instance) -> double { return instance; })
      .def("str", &dolfin::Constant::str);

    //-----------------------------------------------------------------------------
    // dolfin::FacetArea
    py::class_<dolfin::FacetArea, std::shared_ptr<dolfin::FacetArea>,
               dolfin::Expression, dolfin::GenericFunction>
      (m, "FacetArea");

    //-----------------------------------------------------------------------------
    // dolfin::MeshCoordinates
    py::class_<dolfin::MeshCoordinates, std::shared_ptr<dolfin::MeshCoordinates>,
               dolfin::Expression, dolfin::GenericFunction>
      (m, "MeshCoordinates");

    //-----------------------------------------------------------------------------
    // dolfin::Function
    py::class_<dolfin::Function, std::shared_ptr<dolfin::Function>, dolfin::GenericFunction>
      (m, "Function")
      .def(py::init<std::shared_ptr<dolfin::FunctionSpace>>())
      .def(py::init<std::shared_ptr<dolfin::FunctionSpace>, std::shared_ptr<dolfin::GenericVector>>())
      .def("_in", &dolfin::Function::in)
      .def("__call__", [](dolfin::Function& self, std::vector<double>& p)
          {
            // FIXME - remove Array and replace with Eigen in DOLFIN
            const dolfin::Array<double> x(p.size(), p.data());
            Eigen::VectorXd values(self.value_size());
            dolfin::Array<double> _values(self.value_size(), values.data());
            self.eval(_values, x);
            return values;
          })
      .def("extrapolate", &dolfin::Function::extrapolate)
      .def("extrapolate", [](dolfin::Function& instance, const py::object v)
           {
             auto _v = v.attr("_cpp_object").cast<dolfin::Function*>();
             instance.extrapolate(*_v);
           })
      .def("interpolate", (void (dolfin::Function::*)(const dolfin::GenericFunction&))
           &dolfin::Function::interpolate)
      .def("interpolate", [](dolfin::Function& instance, const py::object v)
           {
             auto _v = v.attr("_cpp_object").cast<dolfin::GenericFunction*>();
             instance.interpolate(*_v);
           })
      .def("vector", (std::shared_ptr<dolfin::GenericVector> (dolfin::Function::*)())
           &dolfin::Function::vector);

    // FIXME: why is this floating here?
    m.def("interpolate", [](const dolfin::GenericFunction& f,
                          std::shared_ptr<const dolfin::FunctionSpace> V)
          {
            auto g = std::make_shared<dolfin::Function>(V);
            g->interpolate(f);
            return g;
          });

    // dolfin::FunctionAXPY
    py::class_<dolfin::FunctionAXPY, std::shared_ptr<dolfin::FunctionAXPY>> function_axpy(m, "FunctionAXPY");
    function_axpy
      .def(py::init<std::shared_ptr<const dolfin::Function>, double>())
      .def(py::init<std::shared_ptr<const dolfin::Function>, std::shared_ptr<const dolfin::Function>,
           dolfin::FunctionAXPY::Direction>());

    // dolfin::FunctionAXPY enum
    py::enum_<dolfin::FunctionAXPY::Direction>(function_axpy, "Direction")
      .value("ADD_ADD", dolfin::FunctionAXPY::Direction::ADD_ADD)
      .value("SUB_ADD", dolfin::FunctionAXPY::Direction::SUB_ADD)
      .value("ADD_SUB", dolfin::FunctionAXPY::Direction::ADD_SUB)
      .value("SUB_SUB", dolfin::FunctionAXPY::Direction::SUB_SUB);

    //-----------------------------------------------------------------------------
    // dolfin::FunctionSpace
    py::class_<dolfin::FunctionSpace, std::shared_ptr<dolfin::FunctionSpace>, dolfin::Variable>
      (m, "FunctionSpace", py::dynamic_attr())
      .def(py::init<std::shared_ptr<dolfin::Mesh>, std::shared_ptr<dolfin::FiniteElement>,
           std::shared_ptr<dolfin::GenericDofMap>>())
      .def(py::init<const dolfin::FunctionSpace&>())
      .def("dim", &dolfin::FunctionSpace::dim)
      .def("component", &dolfin::FunctionSpace::component)
      .def("contains", &dolfin::FunctionSpace::contains)
      .def("element", &dolfin::FunctionSpace::element)
      .def("mesh", &dolfin::FunctionSpace::mesh)
      .def("dofmap", &dolfin::FunctionSpace::dofmap)
      .def("set_x", &dolfin::FunctionSpace::set_x)
      .def("sub", (std::shared_ptr<dolfin::FunctionSpace> (dolfin::FunctionSpace::*)(std::size_t) const)
           &dolfin::FunctionSpace::sub)
      .def("tabulate_dof_coordinates", [](const dolfin::FunctionSpace& self)
           {
             const std::size_t gdim = self.element()->geometric_dimension();
             std::vector<double> coords = self.tabulate_dof_coordinates();
             assert(coords.size() % gdim  == 0);

             py::array_t<double> c({coords.size()/gdim, gdim}, coords.data() );
             return c;
           });

  }

}
