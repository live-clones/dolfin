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
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <dolfin/common/Array.h>
#include <dolfin/la/solve.h>
#include <dolfin/la/GenericLinearOperator.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/IndexMap.h>
#include <dolfin/la/LinearAlgebraObject.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/Scalar.h>
#include <dolfin/la/DefaultFactory.h>
#include <dolfin/la/EigenFactory.h>
#include <dolfin/la/EigenMatrix.h>
#include <dolfin/la/EigenVector.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/la/KrylovSolver.h>

#include "../mpi_interface.h"


namespace py = pybind11;


namespace dolfin_wrappers
{
  void la(py::module& m)
  {
    // Note: Do not expose dolfin::LinearAlgebraObject as it has a
    // number of templated member functions
    py::class_<dolfin::LinearAlgebraObject, std::shared_ptr<dolfin::LinearAlgebraObject>, dolfin::Variable>(m, "LinearAlgebraObject");

    // dolfin::IndexMap
    py::class_<dolfin::IndexMap, std::shared_ptr<dolfin::IndexMap>> index_map(m, "IndexMap");
    index_map.def("size", &dolfin::IndexMap::size);

    py::enum_<dolfin::IndexMap::MapSize>(index_map, "MapSize")
      .value("ALL", dolfin::IndexMap::MapSize::ALL)
      .value("OWNED", dolfin::IndexMap::MapSize::OWNED)
      .value("UNOWNED", dolfin::IndexMap::MapSize::UNOWNED)
      .value("GLOBAL", dolfin::IndexMap::MapSize::GLOBAL);

    // dolfin::GenericLinearOperator class
    py::class_<dolfin::GenericLinearOperator, std::shared_ptr<dolfin::GenericLinearOperator>>
      (m, "GenericLinearOperatqor", "DOLFIN GenericLinearOperator object")
      .def("mult", &dolfin::GenericLinearOperator::mult);

    // dolfin::GenericTensor class
    py::class_<dolfin::GenericTensor, std::shared_ptr<dolfin::GenericTensor>, dolfin::LinearAlgebraObject>
      (m, "GenericTensor", "DOLFIN GenericTensor object");

    // dolfin::GenericMatrix class
    py::class_<dolfin::GenericMatrix, std::shared_ptr<dolfin::GenericMatrix>,
               dolfin::GenericTensor, dolfin::GenericLinearOperator>
      (m, "GenericMatrix", "DOLFIN GenericMatrix object")
      .def("init_vector", &dolfin::GenericMatrix::init_vector)
      .def("local_range", &dolfin::GenericMatrix::local_range)
      .def("norm", &dolfin::GenericMatrix::norm)
      .def("nnz", &dolfin::GenericMatrix::nnz)
      .def("size", &dolfin::GenericMatrix::size)
      .def("get_diagonal", &dolfin::GenericMatrix::get_diagonal)
      .def("set_diagonal", &dolfin::GenericMatrix::set_diagonal)
      .def("array", [](const dolfin::GenericMatrix& instance)
           {
             // FIXME: This function is highly dubious. It assumes a
             // particular matrix data layout.

             auto m_range = instance.local_range(0);
             std::size_t num_rows = m_range.second - m_range.first;
             std::size_t num_cols = instance.size(1);

             Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_rows, num_cols);
             std::vector<std::size_t> columns;
             std::vector<double> values;
             for (std::size_t i = 0; i < num_rows; ++i)
             {
               const std::size_t row = i + m_range.first;
               instance.getrow(row, columns, values);
               for (std::size_t j = 0; j < columns.size(); ++j)
                 A(row, columns[j]) = values[j];
             }

             return A;
           });

    // dolfin::GenericVector class
    py::class_<dolfin::GenericVector, std::shared_ptr<dolfin::GenericVector>,
               dolfin::GenericTensor>
      (m, "GenericVector", "DOLFIN GenericVector object")
      .def("__getitem__", [](dolfin::GenericVector& self, py::slice slice)
           {
             std::size_t start, stop, step, slicelength;
             if (!slice.compute(self.size(), &start, &stop, &step, &slicelength))
               throw py::error_already_set();
             if (start != 0 or stop != self.size() or step != 1)
               throw std::range_error("Only full slices are supported");
             std::vector<double> values;
             self.get_local(values);
             return values;
           })
      .def("__getitem__", &dolfin::GenericVector::getitem)
      .def("__setitem__", [](dolfin::GenericVector& self, py::slice slice, double value)
           {
             std::size_t start, stop, step, slicelength;
             if (!slice.compute(self.size(), &start, &stop, &step, &slicelength))
               throw py::error_already_set();
             if (start != 0 or stop != self.size() or step != 1)
               throw std::range_error("Only full slices are supported");
             self = value;
           })
      .def("__setitem__", &dolfin::GenericVector::setitem)
      .def("__sub__", [](dolfin::GenericVector& self, dolfin::GenericVector& v)
           {
             auto a = self.copy();
             (*a) -= v;
             return a;
           })
      .def("get_local", [](const dolfin::GenericVector& instance, const std::vector<long>& rows)
           {
             std::vector<dolfin::la_index> _rows(rows.begin(), rows.end());
             py::array_t<double> data(rows.size());
             instance.get_local(data.mutable_data(), _rows.size(), _rows.data());
             return data;
           })
      .def("set_local", [](dolfin::GenericVector& instance, std::vector<double> values)
           {
             std::vector<dolfin::la_index> indices(values.size());
             std::iota(indices.begin(), indices.end(), 0);
             instance.set_local(values.data(), values.size(), indices.data());
           })
      .def("add_local", [](dolfin::GenericVector& self, py::array_t<double> values)
           {
             std::cout << "Testing shape: " << values.ndim() << ", " << values.size() << std::endl;
             std::cout << values.shape(0) << std::endl;
             assert(values.ndim() == 1);
             dolfin::Array<double> _values(values.size(), values.mutable_data());
             self.add_local(_values);
           })
      .def("gather", [](const dolfin::GenericVector& instance, std::vector<dolfin::la_index> rows)
           {
             std::vector<double> values(rows.size());
             instance.gather(values, rows);
             return py::array_t<double>(values.size(), values.data());
           })
      .def("sum", (double (dolfin::GenericVector::*)() const) &dolfin::GenericVector::sum)
      .def("array", [](const dolfin::GenericVector& instance)
           {
             std::vector<double> values;
             instance.get_local(values);
             return py::array_t<double>(values.size(), values.data());
           });

    // dolfin::Matrix class
    py::class_<dolfin::Matrix, std::shared_ptr<dolfin::Matrix>, dolfin::GenericMatrix>
      (m, "Matrix", "DOLFIN Matrix object", py::multiple_inheritance())
      .def(py::init<>())
      .def(py::init<MPI_Comm>())
      .def("instance", (std::shared_ptr<dolfin::LinearAlgebraObject>(dolfin::Matrix::*)())
           &dolfin::Matrix::shared_instance);

    // dolfin::Vector class
    py::class_<dolfin::Vector, std::shared_ptr<dolfin::Vector>, dolfin::GenericVector>
      (m, "Vector", "DOLFIN Vector object")
      .def(py::init<>())
      .def(py::init<const dolfin::Vector&>())
      .def(py::init<MPI_Comm>())
      .def(py::init<MPI_Comm, std::size_t>())
      .def("init", (void (dolfin::Vector::*)(std::pair<std::size_t, std::size_t>))&dolfin::Vector::init)
      .def("mpi_comm", &dolfin::Vector::mpi_comm)
      .def("size", &dolfin::Vector::size)
      .def("__add__", [](dolfin::Vector& self, dolfin::Vector& v)
           {
             dolfin::Vector a(self);
             a += v;
             return a;
           })
      .def("__add__", [](dolfin::Vector& self, double b)
           {
             dolfin::Vector a(self);
             a += b;
             return a;
           })
      .def("__sub__", [](dolfin::Vector& self, dolfin::Vector& v)
           {
             dolfin::Vector a(self);
             a -= v;
             return a;
           })
      .def("__sub__", [](dolfin::Vector& self, double b)
           {
             dolfin::Vector a(self);
             a -= b;
             return a;
           })
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
      .def("__itruediv__", (const dolfin::Vector& (dolfin::Vector::*)(double))
           &dolfin::Vector::operator/=)
      .def("__setitem__", [](dolfin::Vector& self, dolfin::la_index index, double value)
           { self.instance()->setitem(index, value); })
      .def("__setitem__", [](dolfin::Vector& self, py::slice slice, double value)
           {
             std::size_t start, stop, step, slicelength;
             if (!slice.compute(self.size(), &start, &stop, &step, &slicelength))
               throw py::error_already_set();
             if (start != 0 or stop != self.size() or step != 1)
               throw std::range_error("Only full slices are supported");
             *self.instance() = value; })
      .def("sum", (double (dolfin::Vector::*)() const) &dolfin::Vector::sum)
      .def("min", &dolfin::Vector::min)
      .def("max", &dolfin::Vector::max)
      .def("abs", &dolfin::Vector::abs)
      .def("norm", &dolfin::Vector::norm)
      .def("inner", &dolfin::Vector::inner)
      .def("axpy", &dolfin::Vector::axpy)
      .def("zero", &dolfin::Vector::zero)
      .def("local_size", &dolfin::Vector::local_size)
      .def("local_range", &dolfin::Vector::local_range)
      .def("owns_index", &dolfin::Vector::owns_index)
      .def("get_local", [](dolfin::Vector& self)
           {
             std::vector<double> values;
             self.get_local(values);
             return values;
           })
      .def("apply", &dolfin::Vector::apply)
      .def("str", &dolfin::Vector::str)
      .def("shared_instance", (std::shared_ptr<dolfin::LinearAlgebraObject>(dolfin::Vector::*)())
           &dolfin::Vector::shared_instance);

    //----------------------------------------------------------------------------
    // dolfin::Scalar
    py::class_<dolfin::Scalar, std::shared_ptr<dolfin::Scalar>, dolfin::GenericTensor>
      (m, "Scalar")
      .def(py::init<>())
      .def(py::init<MPI_Comm>());

    //----------------------------------------------------------------------------
    // dolfin::GenericLinearAlgebraFactory class
    py::class_<dolfin::GenericLinearAlgebraFactory, std::shared_ptr<dolfin::GenericLinearAlgebraFactory>>
      (m, "GenericLinearAlgebraFactory", "DOLFIN GenericLinearAlgebraFactory object");

    //----------------------------------------------------------------------------
    // dolfin::DefaultFactory class
    py::class_<dolfin::DefaultFactory, std::shared_ptr<dolfin::DefaultFactory>>
      (m, "DefaultFactory", "DOLFIN DefaultFactory object")
      .def_static("factory", &dolfin::DefaultFactory::factory)
      .def("create_matrix", &dolfin::DefaultFactory::create_matrix)
      .def("create_vector", &dolfin::DefaultFactory::create_vector);

    //----------------------------------------------------------------------------
    // dolfin::EigenFactory class
    py::class_<dolfin::EigenFactory, std::shared_ptr<dolfin::EigenFactory>,
      dolfin::GenericLinearAlgebraFactory>
      (m, "EigenFactory", "DOLFIN EigenFactory object")
      .def("instance", &dolfin::EigenFactory::instance)
      .def("create_matrix", &dolfin::EigenFactory::create_matrix)
      .def("create_vector", &dolfin::EigenFactory::create_vector);

    //----------------------------------------------------------------------------
    // dolfin::EigenVector class
    py::class_<dolfin::EigenVector, std::shared_ptr<dolfin::EigenVector>,
               dolfin::GenericVector>
      (m, "EigenVector", "DOLFIN EigenVector object")
      .def(py::init<MPI_Comm, std::size_t>())
      .def("array", (Eigen::VectorXd& (dolfin::EigenVector::*)()) &dolfin::EigenVector::vec,
           py::return_value_policy::reference_internal);

    //----------------------------------------------------------------------------
    // dolfin::EigenMatrix class
    py::class_<dolfin::EigenMatrix, std::shared_ptr<dolfin::EigenMatrix>,
               dolfin::GenericMatrix>
      (m, "EigenMatrix", "DOLFIN EigenMatrix object")
      .def(py::init<>())
      .def(py::init<std::size_t, std::size_t>())
      .def("array", (dolfin::EigenMatrix::eigen_matrix_type& (dolfin::EigenMatrix::*)()) &dolfin::EigenMatrix::mat,
           py::return_value_policy::reference_internal);

    #ifdef HAS_PETSC
    py::class_<dolfin::PETScObject, std::shared_ptr<dolfin::PETScObject>>(m, "PETScObject");

    // dolfin::PETScFactory class
    py::class_<dolfin::PETScFactory, std::shared_ptr<dolfin::PETScFactory>,
      dolfin::GenericLinearAlgebraFactory>
      (m, "PETScFactory", "DOLFIN PETScFactory object")
      .def("instance", &dolfin::PETScFactory::instance)
      .def("create_matrix", &dolfin::PETScFactory::create_matrix)
      .def("create_vector", &dolfin::PETScFactory::create_vector);

    //----------------------------------------------------------------------------
    // dolfin::PETScVector class
    py::class_<dolfin::PETScVector, std::shared_ptr<dolfin::PETScVector>,
               dolfin::GenericVector, dolfin::PETScObject>
      (m, "PETScVector", "DOLFIN PETScVector object")
      .def(py::init<MPI_Comm>())
      .def(py::init<MPI_Comm, std::size_t>())
      .def("update_ghost_values", &dolfin::PETScVector::update_ghost_values);

    // dolfin::PETScBaseMatrix class
    py::class_<dolfin::PETScBaseMatrix, std::shared_ptr<dolfin::PETScBaseMatrix>,
               dolfin::PETScObject, dolfin::Variable>(m, "PETScBaseMatrix");

    // dolfin::PETScMatrix class
    py::class_<dolfin::PETScMatrix, std::shared_ptr<dolfin::PETScMatrix>,
               dolfin::GenericMatrix, dolfin::PETScBaseMatrix>
      (m, "PETScMatrix", "DOLFIN PETScMatrix object")
      .def(py::init<>())
      .def(py::init<MPI_Comm>());
    #endif

    //-----------------------------------------------------------------------------

    // dolfin::GenericLinearSolver class
    py::class_<dolfin::GenericLinearSolver, std::shared_ptr<dolfin::GenericLinearSolver>,
               dolfin::Variable>
      (m, "GenericLinearSolver", "DOLFIN GenericLinearSolver object");

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

    // Cast to backend type
    m.def("has_type_matrix", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::Matrix>(x); });
    m.def("as_type_matrix", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::Matrix>(x); });

    m.def("has_type_vector", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::Vector>(x); });
    m.def("as_type_vector", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::Vector>(x); });

    m.def("has_type_eigen_matrix", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::EigenMatrix>(x); });
    m.def("as_type_eigen_matrix", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::EigenMatrix>(x); });

    m.def("has_type_eigen_vector", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::EigenVector>(x); });
    m.def("as_type_eigen_vector", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::EigenVector>(x); });
#ifdef HAS_PETSC
    m.def("has_type_petsc_matrix", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::PETScMatrix>(x); });
    m.def("as_type_petsc_matrix", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::PETScMatrix>(x); });

    m.def("has_type_petsc_vector", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::PETScVector>(x); });
    m.def("as_type_petsc_vector", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::PETScVector>(x); });
#endif
#ifdef HAS_TRILINOS
    m.def("has_type_tpetra_matrix", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::TpetraMatrix>(x); });
    m.def("as_type_tpetra_matrix", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::TpetraMatrix>(x); });

    m.def("has_type_tpetra_vector", [](const dolfin::LinearAlgebraObject& x)
          { return dolfin::has_type<dolfin::TpetraVector>(x); });
    m.def("as_type_tpetra_vector", [](std::shared_ptr<dolfin::LinearAlgebraObject> x)
          { return dolfin::as_type<dolfin::TpetraVector>(x); });
#endif

    m.def("has_linear_algebra_backend", &dolfin::has_linear_algebra_backend);
    m.def("linear_algebra_backends", &dolfin::linear_algebra_backends);
    m.def("has_krylov_solver_method", &dolfin::has_krylov_solver_method);
    m.def("has_krylov_solver_preconditioner", &dolfin::has_krylov_solver_preconditioner);

  }
}
