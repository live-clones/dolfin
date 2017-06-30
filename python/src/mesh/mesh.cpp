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

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>

//#include <dolfin/geometry/Point.h>
//#include <dolfin/generation/BoxMesh.h>


namespace py = pybind11;

namespace dolfin_wrappers
{

  void mesh(py::module& m)
  {
    //-----------------------------------------------------------------------------
    // dolfin::Mesh class
    py::class_<dolfin::Mesh, std::shared_ptr<dolfin::Mesh>>
      mesh(m, "Mesh", "DOLFIN Mesh object");

    // Constructors
    // mesh.def(py::init<>());

    mesh.def("num_entities", &dolfin::Mesh::num_entities, "Number of mesh entities");
    mesh.def("topology", (const dolfin::MeshTopology& (dolfin::Mesh::*)() const)
	                                   &dolfin::Mesh::topology, "Mesh topology");
    mesh.def("geometry", (dolfin::MeshGeometry& (dolfin::Mesh::*)())
	                                   &dolfin::Mesh::geometry, "Mesh geometry");
    mesh.def("coordinates",
	     [](const dolfin::Mesh& self)
	     {
	       const std::size_t num_points = self.geometry().num_points();
	       const std::size_t gdim = self.geometry().dim();
	       py::array_t<double, py::array::c_style>
		 f({num_points, gdim}, self.geometry().x().data());
	       return f;
	     });
    
    //-----------------------------------------------------------------------------
    // dolfin::MeshTopology class
    py::class_<dolfin::MeshTopology, std::shared_ptr<dolfin::MeshTopology>>
      mesh_topology(m, "MeshTopology", "DOLFIN MeshTopology object");

    mesh_topology.def("dim", &dolfin::MeshTopology::dim, "Topological dimension");

    //-----------------------------------------------------------------------------
    // dolfin::MeshGeometry class
    py::class_<dolfin::MeshGeometry, std::shared_ptr<dolfin::MeshGeometry>>
      mesh_geometry(m, "MeshGeometry", "DOLFIN MeshGeometry object");

    mesh_geometry.def("dim", &dolfin::MeshGeometry::dim, "Geometrical dimension");
    
    
  }

}
