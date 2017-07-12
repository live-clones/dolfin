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
  // common
  void common(py::module& m);
  void mpi(py::module& m);

  void ale(py::module& m);
  void experimental(py::module& m);
  void fem(py::module& m);
  void function(py::module& m);
  void generation(py::module& m);
  void geometry(py::module& m);
  void graph(py::module& m);
  void io(py::module& m);
  void la(py::module& m);
  void math(py::module& m);
  void mesh(py::module& m);
  void multistage(py::module& m);
  void parameter(py::module& m);
  void refinement(py::module& m);
}


PYBIND11_MODULE(cpp, m)
{
  // Create module for C++ wrappers
  m.doc() ="DOLFIN Python interface";

  // Create MPI class [common]
  dolfin_wrappers::mpi(m);

  // Create ale submodule [ale]
  py::module ale = m.def_submodule("ale", "DOLFIN ALE module");
  dolfin_wrappers::ale(ale);

  // Create common submodule [common]
  py::module common = m.def_submodule("common", "DOLFIN common module");
  dolfin_wrappers::common(common);

  // Create math submodule [math]
  py::module math = m.def_submodule("math", "DOLFIN math library module");
  dolfin_wrappers::math(math);

  // Create mesh submodule [mesh]
  py::module mesh = m.def_submodule("mesh", "DOLFIN mesh library module");
  dolfin_wrappers::mesh(mesh);

  // Create multistage submodule [multistage]
  py::module multistage = m.def_submodule("multistage", "DOLFIN multistage library module");
  dolfin_wrappers::multistage(multistage);

  // Create graph submodule [graph]
  py::module graph = m.def_submodule("graph", "DOLFIN graph module");
  dolfin_wrappers::graph(graph);

  // Create fem submodule [fem]
  py::module fem = m.def_submodule("fem", "DOLFIN FEM module");
  dolfin_wrappers::fem(fem);

  // Create function submodule [function]
  py::module function = m.def_submodule("function",
                                        "DOLFIN function module");
  dolfin_wrappers::function(function);

  // Create generation submodule [generation]
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

  // Create la submodule
  py::module la = m.def_submodule("la", "DOLFIN linear algebra module");
  dolfin_wrappers::la(la);

  // Create parameter submodule
  py::module parameter = m.def_submodule("parameter", "DOLFIN parameter module");
  dolfin_wrappers::parameter(parameter);

  // Create refinement submodule
  py::module refinement = m.def_submodule("refinement", "DOLFIN refinement module");
  dolfin_wrappers::refinement(refinement);

  // Create experimental submodule
  py::module experimental = m.def_submodule("experimental",
                                            "Experimental module");
  dolfin_wrappers::experimental(experimental);
}
