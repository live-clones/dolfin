// Copyright (C) 2017 Chris Richardson
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

#include <dolfin/graph/Graph.h>
#include <dolfin/graph/GraphBuilder.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void graph(py::module& m)
  {
    py::class_<dolfin::Graph>(m, "Graph");

    py::class_<dolfin::GraphBuilder>(m, "GraphBuilder");

  }

}