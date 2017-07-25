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
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dolfin/common/Variable.h>
#include <dolfin/log/log.h>
#include "../mpi_interface.h"

namespace py = pybind11;

namespace dolfin_wrappers
{
  void log(py::module& m)
  {

    //m.def("info", (void (*)(const dolfin::Variable&, bool)) &dolfin::info,
    //      py::arg("variable"), py::arg("verbose")=false);
    m.def("info", [](const dolfin::Variable& v){ dolfin::info(v); });
    m.def("info", [](const dolfin::Variable& v, bool verbose){ dolfin::info(v, verbose); });


  }

}