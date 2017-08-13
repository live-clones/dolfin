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
#include <dolfin/log/Table.h>
#include <dolfin/mesh/Mesh.h>
#include "../mpi_interface.h"

namespace py = pybind11;

namespace dolfin_wrappers
{
  void log(py::module& m)
  {
    py::class_<dolfin::Table, std::shared_ptr<dolfin::Table>>(m, "Table")
      .def(py::init<std::string>())
      .def("str", &dolfin::Table::str);

    //m.def("info", (void (*)(const dolfin::Variable&, bool)) &dolfin::info,
    //      py::arg("variable"), py::arg("verbose")=false);
    m.def("info", [](const dolfin::Variable& v){ dolfin::info(v); });
    m.def("info", [](const dolfin::Variable& v, bool verbose){ dolfin::info(v, verbose); });
    m.def("info", [](std::string s){ dolfin::info(s); });
    m.def("info", [](const dolfin::Parameters& p, bool verbose){ dolfin::info(p, verbose); });
    m.def("info", [](const dolfin::Mesh& mesh, bool verbose){ dolfin::info(mesh, verbose); },
          py::arg("mesh"), py::arg("verbose")=false);

  }

}
