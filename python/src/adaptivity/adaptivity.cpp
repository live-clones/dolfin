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

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dolfin/adaptivity/TimeSeries.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Mesh.h>

#include "../mpi_interface.h"

namespace py = pybind11;

namespace dolfin_wrappers
{
  void adaptivity(py::module& m)
  {
     // Wrap TimesSeries
    py::class_<dolfin::TimeSeries, std::shared_ptr<dolfin::TimeSeries>>(m, "TimeSeries")
      .def(py::init<std::string>())
      .def(py::init<MPI_Comm, std::string>())
      .def("store", (void (dolfin::TimeSeries::*)(const dolfin::GenericVector&, double)) &dolfin::TimeSeries::store)
      .def("store", (void (dolfin::TimeSeries::*)(const dolfin::Mesh&, double)) &dolfin::TimeSeries::store)
      .def("retrieve", (void (dolfin::TimeSeries::*)(dolfin::GenericVector&, double, bool) const) &dolfin::TimeSeries::retrieve,
           py::arg("vector"), py::arg("t"), py::arg("interpolate")=true)
      .def("retrieve", (void (dolfin::TimeSeries::*)(dolfin::Mesh&, double) const) &dolfin::TimeSeries::retrieve)
      .def("vector_times", &dolfin::TimeSeries::vector_times)
      .def("mesh_times", &dolfin::TimeSeries::mesh_times);
  }

}
