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
#include <typeinfo>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <dolfin/la/solve.h>
#include <dolfin/la/GenericLinearOperator.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/LinearAlgebraObject.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/EigenFactory.h>
#include <dolfin/la/EigenMatrix.h>
#include <dolfin/la/EigenVector.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/la/KrylovSolver.h>

#include "../openmpi.h"

namespace py = pybind11;

namespace dolfin_wrappers
{
  void la(py::module& m)
  {
    // dolfin::GenericLinearOperator class
    //py::class_<dolfin::LinearAlgerabObject, std::shared_ptr<dolfin::LinearAlgerabObject>>
    //  (m, "LinearAlgerabObject", "DOLFIN LinearAlgerabObject object");

    // dolfin::GenericLinearOperator class
    py::class_<dolfin::GenericLinearOperator, std::shared_ptr<dolfin::GenericLinearOperator>>
      (m, "GenericLinearOperator", "DOLFIN GenericLinearOperator object");

    // dolfin::GenericTensor class
    py::class_<dolfin::GenericTensor, std::shared_ptr<dolfin::GenericTensor>>
      (m, "GenericTensor", "DOLFIN GenericTensor object");

    // dolfin::GenericMatrix class
    py::class_<dolfin::GenericMatrix, std::shared_ptr<dolfin::GenericMatrix>,
               dolfin::GenericTensor, dolfin::GenericLinearOperator>
      (m, "GenericMatrix", "DOLFIN GenericMatrix object");

    // dolfin::GenericVector class
    py::class_<dolfin::GenericVector, std::shared_ptr<dolfin::GenericVector>,
               dolfin::GenericTensor>
      (m, "GenericVector", "DOLFIN GenericVector object");

    // dolfin::GenericLinearSolver class
    py::class_<dolfin::GenericLinearSolver, std::shared_ptr<dolfin::GenericLinearSolver>>
      (m, "GenericLinearSolver", "DOLFIN GenericLinearSolver object");

    //-----------------------------------------------------------------------------

    // dolfin::Matrix class
    py::class_<dolfin::Matrix, std::shared_ptr<dolfin::Matrix>, dolfin::GenericMatrix,
               dolfin::GenericTensor, dolfin::GenericLinearOperator>
      (m, "Matrix", "DOLFIN Matrix object")
      .def(py::init<>())
      .def(py::init<MPI_Comm>());

    //-----------------------------------------------------------------------------
    // dolfin::Vector class
    py::class_<dolfin::Vector, std::shared_ptr<dolfin::Vector>,
               dolfin::GenericVector, dolfin::GenericTensor>
      (m, "Vector", "DOLFIN Vector object")
      .def(py::init<>())
      .def(py::init<MPI_Comm>())
      .def(py::init<MPI_Comm, std::size_t>())
      .def("__iadd__", (const dolfin::GenericVector& (dolfin::Vector::*)(double))
           &dolfin::Vector::operator+=)
      .def("__iadd__", (const dolfin::Vector& (dolfin::Vector::*)(const dolfin::GenericVector&))
           &dolfin::Vector::operator+=)
      .def("__isub__", (const dolfin::GenericVector& (dolfin::Vector::*)(double))
           &dolfin::Vector::operator-=)
      .def("__isub__", (const dolfin::Vector& (dolfin::Vector::*)(const dolfin::GenericVector&))
           &dolfin::Vector::operator-=)
      .def("__imul__", (const dolfin::Vector& (dolfin::Vector::*)(double))
           &dolfin::Vector::operator*=)
      .def("__imul__", (const dolfin::Vector& (dolfin::Vector::*)(const dolfin::GenericVector&))
           &dolfin::Vector::operator*=)
      .def("__idiv__", (const dolfin::Vector& (dolfin::Vector::*)(double))
           &dolfin::Vector::operator/=)
      .def("__setitem__", [](dolfin::Vector& self, dolfin::la_index index, double value)
           { self.instance()->setitem(index, value); })
      .def("backend_type", [](dolfin::Vector& self)
           {
             // Experiment with determining backend type
             auto instance = self.instance();
             auto type_index = std::type_index(typeid(*instance));
             if (type_index == std::type_index(typeid(dolfin::EigenVector)))
               return "EigenVector";
             else
               return "Not an EigenVector";
           });

    //----------------------------------------------------------------------------
    // dolfin::EigenFactory class
    py::class_<dolfin::EigenFactory, std::shared_ptr<dolfin::EigenFactory>>
      (m, "EigenFactory", "DOLFIN EigenFactory object")
      .def("instance", &dolfin::EigenFactory::instance)
      .def("create_matrix", &dolfin::EigenFactory::create_matrix)
      .def("create_vector", &dolfin::EigenFactory::create_vector);

    //----------------------------------------------------------------------------
    // dolfin::EigenVector class
    py::class_<dolfin::EigenVector, std::shared_ptr<dolfin::EigenVector>,
               dolfin::GenericVector, dolfin::GenericTensor>
      (m, "EigenVector", "DOLFIN EigenVector object")
      .def(py::init<MPI_Comm, std::size_t>())
      .def("array", (Eigen::VectorXd& (dolfin::EigenVector::*)()) &dolfin::EigenVector::vec,
           py::return_value_policy::reference_internal);

    //----------------------------------------------------------------------------
    // dolfin::EigenMatrix class
    py::class_<dolfin::EigenMatrix, std::shared_ptr<dolfin::EigenMatrix>,
               dolfin::GenericMatrix, dolfin::GenericTensor, dolfin::GenericLinearOperator>
      (m, "EigenMatrix", "DOLFIN EigenMatrix object")
      .def(py::init<>())
      .def(py::init<std::size_t, std::size_t>())
      .def("array", (dolfin::EigenMatrix::eigen_matrix_type& (dolfin::EigenMatrix::*)()) &dolfin::EigenMatrix::mat,
           py::return_value_policy::reference_internal);

    #ifdef HAS_PETSC
    //----------------------------------------------------------------------------
    // dolfin::PETScVector class
    py::class_<dolfin::PETScVector, std::shared_ptr<dolfin::PETScVector>,
               dolfin::GenericVector, dolfin::GenericTensor>
      (m, "PETScVector", "DOLFIN PETScVector object")
      .def(py::init<MPI_Comm>())
      .def(py::init<MPI_Comm, std::size_t>());

    //----------------------------------------------------------------------------
    // dolfin::PETScMatrix class
    py::class_<dolfin::PETScMatrix, std::shared_ptr<dolfin::PETScMatrix>,
               dolfin::GenericMatrix, dolfin::GenericTensor,
               dolfin::GenericLinearOperator>
      (m, "PETScMatrix", "DOLFIN PETScMatrix object")
      .def(py::init<>())
      .def(py::init<MPI_Comm>());
    #endif

    //-----------------------------------------------------------------------------
    // dolfin::LUSolver class
    py::class_<dolfin::LUSolver, std::shared_ptr<dolfin::LUSolver>>
      (m, "LUSolver", "DOLFIN LUSolver object")
      .def(py::init<MPI_Comm, std::shared_ptr<const dolfin::GenericLinearOperator>,
           std::string>(),
           py::arg("comm"), py::arg("A"), py::arg("method") = "default")
      .def("solve", (std::size_t (dolfin::LUSolver::*)(dolfin::GenericVector&,
                                                       const dolfin::GenericVector&))
           &dolfin::LUSolver::solve);

    //-----------------------------------------------------------------------------
    // dolfin::KrylovSolver class
    py::class_<dolfin::KrylovSolver, std::shared_ptr<dolfin::KrylovSolver>,
               dolfin::GenericLinearSolver>
      (m, "KrylovSolver", "DOLFIN KrylovSolver object")
      .def(py::init<std::shared_ptr<const dolfin::GenericLinearOperator>,
           std::string, std::string>(), py::arg("A"),
           py::arg("method")="default", py::arg("preconditioner")="default")
      .def(py::init<MPI_Comm, std::shared_ptr<const dolfin::GenericLinearOperator>,
           std::string, std::string>(), py::arg("comm"), py::arg("A"),
           py::arg("method")="default", py::arg("preconditioner")="default")
      .def("solve", (std::size_t (dolfin::KrylovSolver::*)(dolfin::GenericVector&,
                                                           const dolfin::GenericVector&))
           &dolfin::KrylovSolver::solve);


    m.def("has_linear_algebra_backend", &dolfin::has_linear_algebra_backend);
    m.def("linear_algebra_backends", &dolfin::linear_algebra_backends);
    m.def("has_krylov_solver_method", &dolfin::has_krylov_solver_method);
    m.def("has_krylov_solver_preconditioner", &dolfin::has_krylov_solver_preconditioner);
  }
}