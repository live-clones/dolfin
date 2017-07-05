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

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <dolfin/common/MPI.h>
#include <dolfin/common/Variable.h>

namespace py = pybind11;

namespace dolfin_wrappers
{
  void common(py::module& m)
  {
    // Variable
    py::class_<dolfin::Variable, std::shared_ptr<dolfin::Variable>>
      (m, "Variable", "Variable base class")
      .def("id", &dolfin::Variable::id)
      .def("name", &dolfin::Variable::name)
      .def("rename", &dolfin::Variable::rename);
  }

  void mpi(py::module& m)
  {
    // MPI
    m.attr("comm_world") = MPI_COMM_WORLD;
    m.attr("comm_self") = MPI_COMM_SELF;

    m.def("rank", &dolfin::MPI::rank);
    m.def("size", &dolfin::MPI::size);
    m.def("max", &dolfin::MPI::max<double>);
    m.def("min", &dolfin::MPI::min<double>);
    m.def("sum", &dolfin::MPI::sum<double>);

  }


}
