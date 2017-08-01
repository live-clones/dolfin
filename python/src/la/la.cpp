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
#include <pybind11/operators.h>

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
#include <dolfin/la/TensorLayout.h>
#include <dolfin/la/DefaultFactory.h>
#include <dolfin/la/EigenFactory.h>
#include <dolfin/la/EigenMatrix.h>
#include <dolfin/la/EigenVector.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScOptions.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/la/KrylovSolver.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/solve.h>

#include "../mpi_interface.h"


namespace py = pybind11;


namespace dolfin_wrappers
{
  void la(py::module& m)
  {
    // dolfin::IndexMap
    py::class_<dolfin::IndexMap, std::shared_ptr<dolfin::IndexMap>> index_map(m, "IndexMap");
    index_map.def("size", &dolfin::IndexMap::size);

    // dolfin::IndexMap enums
    py::enum_<dolfin::IndexMap::MapSize>(index_map, "MapSize")
      .value("ALL", dolfin::IndexMap::MapSize::ALL)
      .value("OWNED", dolfin::IndexMap::MapSize::OWNED)
      .value("UNOWNED", dolfin::IndexMap::MapSize::UNOWNED)
      .value("GLOBAL", dolfin::IndexMap::MapSize::GLOBAL);

    // dolfin::SparsityPattern
    py::class_<dolfin::SparsityPattern, std::shared_ptr<dolfin::SparsityPattern>>(m, "SparsityPattern")
      .def("init", &dolfin::SparsityPattern::init)
      .def("num_nonzeros", &dolfin::SparsityPattern::num_nonzeros)
      .def("num_nonzeros_diagonal", [](const dolfin::SparsityPattern& instance)
           {
             std::vector<std::size_t> num_nonzeros;
             instance.num_nonzeros_diagonal(num_nonzeros);
             return py::array_t<std::size_t>(num_nonzeros.size(), num_nonzeros.data());
           })
      .def("num_nonzeros_off_diagonal", [](const dolfin::SparsityPattern& instance)
           {
             std::vector<std::size_t> num_nonzeros;
             instance.num_nonzeros_off_diagonal(num_nonzeros);
             return py::array_t<std::size_t>(num_nonzeros.size(), num_nonzeros.data());
           })
      .def("num_local_nonzeros", [](const dolfin::SparsityPattern& instance)
           {
             std::vector<std::size_t> num_nonzeros;
             instance.num_local_nonzeros(num_nonzeros);
             return py::array_t<std::size_t>(num_nonzeros.size(), num_nonzeros.data());
           });

    // dolfin::TensorLayout
    py::class_<dolfin::TensorLayout, std::shared_ptr<dolfin::TensorLayout>> tensor_layout(m, "TensorLayout");

    // dolfin::TensorLayout enums
    py::enum_<dolfin::TensorLayout::Sparsity>(tensor_layout, "Sparsity")
      .value("SPARSE", dolfin::TensorLayout::Sparsity::SPARSE)
      .value("DENSE", dolfin::TensorLayout::Sparsity::DENSE);
    py::enum_<dolfin::TensorLayout::Ghosts>(tensor_layout, "Ghosts")
      .value("GHOSTED", dolfin::TensorLayout::Ghosts::GHOSTED)
      .value("UNGHOSTED", dolfin::TensorLayout::Ghosts::UNGHOSTED);

    tensor_layout
      .def(py::init<MPI_Comm, std::size_t, dolfin::TensorLayout::Sparsity>())
      .def(py::init<MPI_Comm, std::vector<std::shared_ptr<const dolfin::IndexMap>>,
           std::size_t, dolfin::TensorLayout::Sparsity, dolfin::TensorLayout::Ghosts>())
      .def("init", &dolfin::TensorLayout::init)
      .def("sparsity_pattern", (std::shared_ptr<dolfin::SparsityPattern> (dolfin::TensorLayout::*)()) &dolfin::TensorLayout::sparsity_pattern);

    // dolfin::LinearAlgebraObject
    py::class_<dolfin::LinearAlgebraObject, std::shared_ptr<dolfin::LinearAlgebraObject>,
               dolfin::Variable>(m, "LinearAlgebraObject");

    // dolfin::GenericLinearOperator class
    py::class_<dolfin::GenericLinearOperator, std::shared_ptr<dolfin::GenericLinearOperator>,
               dolfin::LinearAlgebraObject>
      (m, "GenericLinearOperator", "DOLFIN GenericLinearOperator object")
      .def("mult", &dolfin::GenericLinearOperator::mult);

    // dolfin::GenericTensor class
    py::class_<dolfin::GenericTensor, std::shared_ptr<dolfin::GenericTensor>,
               dolfin::LinearAlgebraObject>
      (m, "GenericTensor", "DOLFIN GenericTensor object")
      .def("init", &dolfin::GenericTensor::init)
      .def("zero", &dolfin::GenericTensor::zero);

    // dolfin::GenericMatrix class
    py::class_<dolfin::GenericMatrix, std::shared_ptr<dolfin::GenericMatrix>,
               dolfin::GenericTensor, dolfin::GenericLinearOperator>
      (m, "GenericMatrix", "DOLFIN GenericMatrix object")
      .def("init_vector", &dolfin::GenericMatrix::init_vector)
      .def("transpmult", &dolfin::GenericMatrix::transpmult)
      .def("__mul__", [](const dolfin::GenericMatrix& self, const dolfin::GenericVector& x)
           {
             dolfin::Vector y;
             self.init_vector(y, 0);
             self.mult(x, y);
             return y;
           }, py::is_operator())
      .def("copy", &dolfin::GenericMatrix::copy)
      .def("local_range", &dolfin::GenericMatrix::local_range)
      .def("norm", &dolfin::GenericMatrix::norm)
      .def("nnz", &dolfin::GenericMatrix::nnz)
      .def("size", &dolfin::GenericMatrix::size)
      .def("get_diagonal", &dolfin::GenericMatrix::get_diagonal)
      .def("set_diagonal", &dolfin::GenericMatrix::set_diagonal)
      .def("getrow", [](const dolfin::GenericMatrix& instance, std::size_t row)
           {
             std::vector<double> values;
             std::vector<std::size_t> columns;
             instance.getrow(row, columns, values);
             return std::make_pair(columns, values);
           })
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
      .def("init", (void (dolfin::GenericVector::*)(std::size_t)) &dolfin::GenericVector::init)
      .def("init", (void (dolfin::GenericVector::*)(const dolfin::TensorLayout&)) &dolfin::GenericVector::init)
      .def("init", (void (dolfin::GenericVector::*)(std::pair<std::size_t, std::size_t>)) &dolfin::GenericVector::init)
      .def("copy", &dolfin::GenericVector::copy)
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
           }, py::is_operator())
      .def("__len__", [](dolfin::GenericVector& self) { return self.size(); })
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
      .def("norm", &dolfin::GenericVector::norm)
      .def("array", [](const dolfin::GenericVector& instance)
           {
             std::vector<double> values;
             instance.get_local(values);
             return py::array_t<double>(values.size(), values.data());
           });

    // dolfin::Matrix class
    py::class_<dolfin::Matrix, std::shared_ptr<dolfin::Matrix>, dolfin::GenericMatrix>
      (m, "Matrix", "DOLFIN Matrix object")
      .def(py::init<>())
      .def(py::init<const dolfin::Matrix&>())  // Remove? (use copy instead)
      .def(py::init<const dolfin::GenericMatrix&>())  // Remove? (use copy instead)
      .def(py::init<MPI_Comm>()) // This comes last of constructors so pybind11 attempts it lasts (avoid OpenMPI comm casting problems)
      .def("instance", (std::shared_ptr<dolfin::LinearAlgebraObject>(dolfin::Matrix::*)())
           &dolfin::Matrix::shared_instance);



    // dolfin::Vector class
    py::class_<dolfin::Vector, std::shared_ptr<dolfin::Vector>, dolfin::GenericVector>
      (m, "Vector", "DOLFIN Vector object")
      .def(py::init<>())
      .def(py::init<const dolfin::Vector&>())
      .def(py::init<MPI_Comm>())
      .def(py::init<MPI_Comm, std::size_t>())
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
      .def(py::self += py::self)
      //      .def("__iadd__", (const dolfin::Vector& (dolfin::Vector::*)(const dolfin::GenericVector&))
      //           &dolfin::Vector::operator+=)
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
      .def("instance", (std::shared_ptr<dolfin::LinearAlgebraObject>(dolfin::Vector::*)())
           &dolfin::Vector::shared_instance);

    //--------------------------------------------------------------------------
    // dolfin::Scalar
    py::class_<dolfin::Scalar, std::shared_ptr<dolfin::Scalar>, dolfin::GenericTensor>
      (m, "Scalar")
      .def(py::init<>())
      .def(py::init<MPI_Comm>())
      .def("add_local_value", &dolfin::Scalar::add_local_value)
      .def("apply", &dolfin::Scalar::apply)
      .def("mpi_comm", &dolfin::Scalar::mpi_comm)
      .def("get_scalar_value", &dolfin::Scalar::get_scalar_value);

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
      .def("sparray", (dolfin::EigenMatrix::eigen_matrix_type& (dolfin::EigenMatrix::*)()) &dolfin::EigenMatrix::mat,
           py::return_value_policy::reference_internal)
      .def("data_view", [](dolfin::EigenMatrix& instance)
           {
             auto _data = instance.data();
             std::size_t nnz = std::get<3>(_data);

             Eigen::Map<const Eigen::VectorXi> rows(std::get<0>(_data), instance.size(0) + 1);
             Eigen::Map<const Eigen::VectorXi> cols(std::get<1>(_data), nnz);
             Eigen::Map<const Eigen::VectorXd> values(std::get<2>(_data), nnz);

             return py::make_tuple(rows, cols, values);
           },
           py::return_value_policy::reference_internal, "Return CSR matrix data as NumPy arrays (shared data)")
      .def("data", [](dolfin::EigenMatrix& instance)
           {
             auto _data = instance.data();
             std::size_t nnz = std::get<3>(_data);

             Eigen::VectorXi rows = Eigen::Map<const Eigen::VectorXi>(std::get<0>(_data), instance.size(0) + 1);
             Eigen::VectorXi cols = Eigen::Map<const Eigen::VectorXi>(std::get<1>(_data), nnz);
             Eigen::VectorXd values  = Eigen::Map<const Eigen::VectorXd>(std::get<2>(_data), nnz);

             return py::make_tuple(rows, cols, values);
           },
           py::return_value_policy::copy, "Return copy of CSR matrix data as NumPy arrays");

    #ifdef HAS_PETSC
    py::class_<dolfin::PETScOptions>(m, "PETScOptions")
      .def_static("set", (void (*)(std::string)) &dolfin::PETScOptions::set)
      .def_static("set", (void (*)(std::string, bool)) &dolfin::PETScOptions::set)
      .def_static("set", (void (*)(std::string, int)) &dolfin::PETScOptions::set)
      .def_static("set", (void (*)(std::string, double)) &dolfin::PETScOptions::set)
      .def_static("set", (void (*)(std::string, std::string)) &dolfin::PETScOptions::set)
      .def_static("clear", (void (*)(std::string)) &dolfin::PETScOptions::clear)
      .def_static("clear", (void (*)()) &dolfin::PETScOptions::clear);

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

    m.def("has_linear_algebra_backend", &dolfin::has_linear_algebra_backend);
    m.def("linear_algebra_backends", &dolfin::linear_algebra_backends);
    m.def("has_krylov_solver_method", &dolfin::has_krylov_solver_method);
    m.def("has_krylov_solver_preconditioner", &dolfin::has_krylov_solver_preconditioner);

    // normalize
    m.def("normalize", &dolfin::normalize);

    // solve
    m.def("solve", (std::size_t (*)(const dolfin::GenericLinearOperator&, dolfin::GenericVector&,
                                    const dolfin::GenericVector&, std::string, std::string)) &dolfin::solve,
          py::arg("A"), py::arg("x"), py::arg("b"), py::arg("method")="lu",
          py::arg("preconditioner")="none");

  }
}
