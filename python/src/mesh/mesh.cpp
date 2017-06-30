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

#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void mesh(py::module& m)
  {
    //-----------------------------------------------------------------------------
    // dolfin::Mesh class
    py::class_<dolfin::Mesh, std::shared_ptr<dolfin::Mesh>>(m, "Mesh", py::dynamic_attr(), "DOLFIN Mesh object")
      .def("num_entities", &dolfin::Mesh::num_entities, "Number of mesh entities")
      .def("topology", (const dolfin::MeshTopology& (dolfin::Mesh::*)() const)
           &dolfin::Mesh::topology, "Mesh topology")
      .def("geometry", (dolfin::MeshGeometry& (dolfin::Mesh::*)())
           &dolfin::Mesh::geometry, "Mesh geometry")
      .def("coordinates",
           [](dolfin::Mesh& self)
           { return py::array({self.geometry().num_points(), self.geometry().dim()},
                              self.geometry().x().data()); })
      .def("coordinates_vec",
           [](dolfin::Mesh& self){ return self.geometry().x(); },
           py::return_value_policy::reference)
      .def("cells",
           [](const dolfin::Mesh& self)
           {
             const unsigned int tdim = self.topology().dim();
             return py::array({self.topology().size(tdim),
                   self.type().num_vertices(tdim)},
               self.topology()(tdim, 0)().data());
           })
      // UFL related
      .def("ufl_id", [](const dolfin::Mesh& self){ return self.id(); })
      .def("cell_name", [](const dolfin::Mesh& self)
           {
             auto gdim = self.geometry().dim();
             auto cellname = self.type().description(false);
             return cellname;
           }
        );

    mesh.def("cell_type",
             [](dolfin::Mesh& self)
             {
               return dolfin::CellType::type2string(self.type().cell_type());
             });

    //-----------------------------------------------------------------------------
    // dolfin::MeshTopology class
    py::class_<dolfin::MeshTopology, std::shared_ptr<dolfin::MeshTopology>>
      mesh_topology(m, "MeshTopology", "DOLFIN MeshTopology object");

    mesh_topology.def("dim", &dolfin::MeshTopology::dim, "Topological dimension");

    //-----------------------------------------------------------------------------
    // dolfin::MeshGeometry class
    py::class_<dolfin::MeshGeometry, std::shared_ptr<dolfin::MeshGeometry>>
      (m, "MeshGeometry", "DOLFIN MeshGeometry object")
      .def("dim", &dolfin::MeshGeometry::dim, "Geometrical dimension")
      .def("degree", &dolfin::MeshGeometry::degree, "Degree");


  }

}
