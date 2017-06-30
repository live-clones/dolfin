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


#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace dolfin_wrappers
{
  void mesh(py::module& m);
  void experimental(py::module& m);
  void fem(py::module& m);
  void function(py::module& m);
  void generation(py::module& m);
  void geometry(py::module& m);
  void io(py::module& m);
}


PYBIND11_MODULE(cpp, m)
{
  // Create module
  m.doc() ="DOLFIN Python interface";

  // Create mesh submodule
  py::module mesh = m.def_submodule("mesh", "DOLFIN mesh library");
  dolfin_wrappers::mesh(mesh);

  // Create fem submodule
  py::module fem = m.def_submodule("fem", "DOLFIN FEM module");
  dolfin_wrappers::fem(fem);

  // Create function submodule
  py::module function = m.def_submodule("function",
                                        "DOLFIN function module");
  dolfin_wrappers::function(function);

  // Create generation submodule
  py::module generation = m.def_submodule("generation",
                                          "DOLFIN mesh generation module");
  dolfin_wrappers::generation(generation);

  // Create geometry submodule
  py::module geometry = m.def_submodule("geometry",
                                        "DOLFIN geometry module");
  dolfin_wrappers::geometry(geometry);

  // Create io submodule
  py::module io = m.def_submodule("io", "DOLFIN I/O module");
  dolfin_wrappers::io(io);

  // Create experimental submodule
  py::module experimental = m.def_submodule("experimental",
                                            "Experimental module");
  dolfin_wrappers::experimental(experimental);
}
