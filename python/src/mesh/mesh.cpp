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
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/SubDomain.h>

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
      .def("coordinates", [](dolfin::Mesh& self)
           {
             return Eigen::Map<Eigen::MatrixXd>(self.geometry().x().data(),
                                                self.geometry().num_points(),
                                                self.geometry().dim());
           })
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
             return dolfin::CellType::type2string(self.type().cell_type());
           }
        );

    //-------------------------------------------------------------------------
    // dolfin::MeshTopology class
    py::class_<dolfin::MeshTopology, std::shared_ptr<dolfin::MeshTopology>>
      (m, "MeshTopology", "DOLFIN MeshTopology object")
      .def("dim", &dolfin::MeshTopology::dim, "Topological dimension");

    //--------------------------------------------------------------------------
    // dolfin::MeshGeometry class
    py::class_<dolfin::MeshGeometry, std::shared_ptr<dolfin::MeshGeometry>>
      (m, "MeshGeometry", "DOLFIN MeshGeometry object")
      .def("dim", &dolfin::MeshGeometry::dim, "Geometrical dimension")
      .def("degree", &dolfin::MeshGeometry::degree, "Degree");

    //--------------------------------------------------------------------------
    // dolfin::SubDomain class

    // FIXME: move somewhere else
    py::class_<dolfin::Array<double>, std::shared_ptr<dolfin::Array<double>>> (m, "Array")
      .def("__init__", [](dolfin::Array<double>& instance, std::vector<double>& x)
           {
             new (&instance) dolfin::Array<double>(x.size(), x.data());
           })
      .def("__getitem__", [](dolfin::Array<double>& self, int i)
           { return self[i]; });

    class PySubDomain : public dolfin::SubDomain
    {
      using dolfin::SubDomain::SubDomain;

      bool inside(const Eigen::Ref<Eigen::VectorXd>& x, bool on_boundary) const override
      {
        PYBIND11_OVERLOAD(bool, dolfin::SubDomain, inside, x, on_boundary);
      }
    };

    py::class_<dolfin::SubDomain, std::shared_ptr<dolfin::SubDomain>, PySubDomain>
      (m, "SubDomain", "DOLFIN SubDomain object")
      .def(py::init<>())
      .def("inside", (bool (dolfin::SubDomain::*)(const Eigen::Ref<Eigen::VectorXd>&, bool) const)
           &dolfin::SubDomain::inside)
      .def("test", [](dolfin::SubDomain& self, bool sw, int n)
           {
             std::vector<double> data(n*3);
             if(sw)
             {
               for (unsigned int i = 0; i < n; ++i)
               {
                 Eigen::Map<Eigen::VectorXd> x(&data[i*3], 3);
               }
             }
             else
             {
               for (unsigned int i = 0; i < n; ++i)
               {
                 dolfin::Array<double> x(3, &data[i*3]);
               }
             }
             return 0;
           });

  }

}
