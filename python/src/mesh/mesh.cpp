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
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>

#include "../openmpi.h"

namespace py = pybind11;

namespace dolfin_wrappers
{

  void mesh(py::module& m)
  {
    m.def("make_dolfin_subdomain",
          [](std::uintptr_t e)
          {
            dolfin::SubDomain *p = reinterpret_cast<dolfin::SubDomain *>(e);
            return std::shared_ptr<const dolfin::SubDomain>(p);
          });

    //-----------------------------------------------------------------------------
    // dolfin::Mesh class
    py::class_<dolfin::Mesh, std::shared_ptr<dolfin::Mesh>>(m, "Mesh", py::dynamic_attr(), "DOLFIN Mesh object")
      .def(py::init<>())
      .def("cells",
           [](const dolfin::Mesh& self)
           {
             const unsigned int tdim = self.topology().dim();
             return py::array({self.topology().size(tdim),
                   self.type().num_vertices(tdim)},
               self.topology()(tdim, 0)().data());
           })
      .def("cell_orientations", &dolfin::Mesh::cell_orientations)
      .def("coordinates", [](dolfin::Mesh& self)
           {
             return Eigen::Map<Eigen::MatrixXd>(self.geometry().x().data(),
                                                self.geometry().num_points(),
                                                self.geometry().dim());
           })
      .def("geometry", (dolfin::MeshGeometry& (dolfin::Mesh::*)())
           &dolfin::Mesh::geometry, "Mesh geometry")
      .def("init_global", &dolfin::Mesh::init_global)
      .def("init", (void (dolfin::Mesh::*)() const) &dolfin::Mesh::init)
      .def("init", (std::size_t (dolfin::Mesh::*)(std::size_t) const) &dolfin::Mesh::init)
      .def("init", (void (dolfin::Mesh::*)(std::size_t, std::size_t) const) &dolfin::Mesh::init)
      .def("init_cell_orientations", &dolfin::Mesh::init_cell_orientations)
      .def("mpi_comm", &dolfin::Mesh::mpi_comm)
      .def("num_entities", &dolfin::Mesh::num_entities, "Number of mesh entities")
      .def("num_vertices", &dolfin::Mesh::num_vertices, "Number of vertices")
      .def("num_cells", &dolfin::Mesh::num_cells, "Number of cells")
      .def("rmin", &dolfin::Mesh::rmin)
      .def("size_global", &dolfin::Mesh::size_global)
      .def("topology", (const dolfin::MeshTopology& (dolfin::Mesh::*)() const)
           &dolfin::Mesh::topology, "Mesh topology")
      // UFL related
      .def("ufl_id", [](const dolfin::Mesh& self){ return self.id(); })
      .def("cell_name", [](const dolfin::Mesh& self)
           { return dolfin::CellType::type2string(self.type().cell_type()); }
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

    //-------------------------------------------------------------------------
    // dolfin::MeshEntity class
    py::class_<dolfin::MeshEntity, std::shared_ptr<dolfin::MeshEntity>>
      (m, "MeshEntity", "DOLFIN MeshEntity object")
      .def(py::init<const dolfin::Mesh&, std::size_t, std::size_t>())
      .def("mesh", &dolfin::MeshEntity::mesh, "Associated mesh")
      .def("dim", &dolfin::MeshEntity::dim, "Topological dimension")
      .def("index", (std::size_t (dolfin::MeshEntity::*)() const)
           &dolfin::MeshEntity::index, "Index")
      .def("global_index", &dolfin::MeshEntity::global_index, "Global index")
      .def("num_entities", &dolfin::MeshEntity::num_entities,
           "Number of incident entities of given dimension")
      .def("num_global_entities", &dolfin::MeshEntity::num_global_entities,
           "Global number of incident entities of given dimension")
      .def("entities", [](dolfin::MeshEntity& self, std::size_t dim)
        {
          return Eigen::Map<const Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>>
            (self.entities(dim), self.num_entities(dim));
        })
      .def("midpoint", &dolfin::MeshEntity::midpoint, "Midpoint of Entity")
      .def("__str__", [](dolfin::MeshEntity& self){return self.str(false);});

    //--------------------------------------------------------------------------
    // dolfin::Vertex class
    py::class_<dolfin::Vertex, std::shared_ptr<dolfin::Vertex>, dolfin::MeshEntity>
      (m, "Vertex", "DOLFIN Vertex object")
      .def(py::init<const dolfin::Mesh&, std::size_t>());

    //--------------------------------------------------------------------------
    // dolfin::Edge class
    py::class_<dolfin::Edge, std::shared_ptr<dolfin::Edge>, dolfin::MeshEntity>
      (m, "Edge", "DOLFIN Edge object")
      .def(py::init<const dolfin::Mesh&, std::size_t>())
      .def("dot", &dolfin::Edge::dot)
      .def("length", &dolfin::Edge::length);

    //--------------------------------------------------------------------------
    // dolfin::Face class
    py::class_<dolfin::Face, std::shared_ptr<dolfin::Face>, dolfin::MeshEntity>
      (m, "Face", "DOLFIN Face object")
      .def(py::init<const dolfin::Mesh&, std::size_t>())
      .def("normal", (dolfin::Point (dolfin::Face::*)() const) &dolfin::Face::normal)
      .def("normal", (double (dolfin::Face::*)(std::size_t) const) &dolfin::Face::normal)
      .def("area", &dolfin::Face::area);

    //--------------------------------------------------------------------------
    // dolfin::Facet class
    py::class_<dolfin::Facet, std::shared_ptr<dolfin::Facet>, dolfin::MeshEntity>
      (m, "Facet", "DOLFIN Facet object")
      .def(py::init<const dolfin::Mesh&, std::size_t>());

    //--------------------------------------------------------------------------
    // dolfin::Cell class
    py::class_<dolfin::Cell, std::shared_ptr<dolfin::Cell>, dolfin::MeshEntity>
      (m, "Cell", "DOLFIN Cell object")
      .def(py::init<const dolfin::Mesh&, std::size_t>())
      .def("distance", &dolfin::Cell::distance)
      .def("facet_area", &dolfin::Cell::facet_area)
      .def("h", &dolfin::Cell::h)
      .def("inradius", &dolfin::Cell::inradius)
      .def("circumradius", &dolfin::Cell::circumradius)
      .def("radius_ratio", &dolfin::Cell::radius_ratio)
      .def("volume", &dolfin::Cell::volume);

    //--------------------------------------------------------------------------
    // dolfin::MeshEntityIterator class
    py::class_<dolfin::MeshEntityIterator, std::shared_ptr<dolfin::MeshEntityIterator>>
      (m, "MeshEntityIterator", "DOLFIN MeshEntityIterator object")
      .def(py::init<const dolfin::Mesh&, std::size_t>())
      .def("__iter__",[](dolfin::MeshEntityIterator& self) { self.operator--(); return self; })
      .def("__next__",[](dolfin::MeshEntityIterator& self) {
          self.operator++();
          if (self.end())
            throw py::stop_iteration("");
          return *self;
        });


    // FIXME: avoid repeating code (macro?)

    py::class_<dolfin::MeshEntityIteratorBase<dolfin::Cell>,
               std::shared_ptr<dolfin::MeshEntityIteratorBase<dolfin::Cell>>>
      (m, "CellIterator", "DOLFIN CellIterator object")
      .def(py::init<const dolfin::Mesh&>())
      .def("__iter__",[](dolfin::MeshEntityIteratorBase<dolfin::Cell>& self) { self.operator--(); return self; })
      .def("__next__",[](dolfin::MeshEntityIteratorBase<dolfin::Cell>& self) {
          self.operator++();
          if (self.end())
            throw py::stop_iteration("");
          return *self;
        });

    m.def("cells", [](dolfin::Mesh& mesh)
          { return dolfin::MeshEntityIteratorBase<dolfin::Cell>(mesh); });
    m.def("cells", [](dolfin::MeshEntity& meshentity)
          { return dolfin::MeshEntityIteratorBase<dolfin::Cell>(meshentity); });

    py::class_<dolfin::MeshEntityIteratorBase<dolfin::Facet>,
               std::shared_ptr<dolfin::MeshEntityIteratorBase<dolfin::Facet>>>
      (m, "FacetIterator", "DOLFIN FacetIterator object")
      .def(py::init<const dolfin::Mesh&>())
      .def("__iter__",[](dolfin::MeshEntityIteratorBase<dolfin::Facet>& self) { self.operator--(); return self; })
      .def("__next__",[](dolfin::MeshEntityIteratorBase<dolfin::Facet>& self) {
          self.operator++();
          if (self.end())
            throw py::stop_iteration("");
          return *self;
        });

    m.def("facets", [](dolfin::Mesh& mesh)
          { return dolfin::MeshEntityIteratorBase<dolfin::Facet>(mesh); });
    m.def("facets", [](dolfin::MeshEntity& meshentity)
          { return dolfin::MeshEntityIteratorBase<dolfin::Facet>(meshentity); });

    py::class_<dolfin::MeshEntityIteratorBase<dolfin::Face>,
               std::shared_ptr<dolfin::MeshEntityIteratorBase<dolfin::Face>>>
      (m, "FaceIterator", "DOLFIN FaceIterator object")
      .def(py::init<const dolfin::Mesh&>())
      .def("__iter__",[](dolfin::MeshEntityIteratorBase<dolfin::Face>& self) { self.operator--(); return self; })
      .def("__next__",[](dolfin::MeshEntityIteratorBase<dolfin::Face>& self) {
          self.operator++();
          if (self.end())
            throw py::stop_iteration("");
          return *self;
        });

    m.def("faces", [](dolfin::Mesh& mesh)
          { return dolfin::MeshEntityIteratorBase<dolfin::Face>(mesh); });
    m.def("faces", [](dolfin::MeshEntity& meshentity)
          { return dolfin::MeshEntityIteratorBase<dolfin::Face>(meshentity); });

    py::class_<dolfin::MeshEntityIteratorBase<dolfin::Edge>,
               std::shared_ptr<dolfin::MeshEntityIteratorBase<dolfin::Edge>>>
      (m, "EdgeIterator", "DOLFIN EdgeIterator object")
      .def(py::init<const dolfin::Mesh&>())
      .def("__iter__",[](dolfin::MeshEntityIteratorBase<dolfin::Edge>& self) { self.operator--(); return self; })
      .def("__next__",[](dolfin::MeshEntityIteratorBase<dolfin::Edge>& self) {
          self.operator++();
          if (self.end())
            throw py::stop_iteration("");
          return *self;
        });

    m.def("edges", [](dolfin::Mesh& mesh)
          { return dolfin::MeshEntityIteratorBase<dolfin::Edge>(mesh); });
    m.def("edges", [](dolfin::MeshEntity& meshentity)
          { return dolfin::MeshEntityIteratorBase<dolfin::Edge>(meshentity); });

    py::class_<dolfin::MeshEntityIteratorBase<dolfin::Vertex>,
               std::shared_ptr<dolfin::MeshEntityIteratorBase<dolfin::Vertex>>>
      (m, "VertexIterator", "DOLFIN VertexIterator object")
      .def(py::init<const dolfin::Mesh&>())
      .def("__iter__",[](dolfin::MeshEntityIteratorBase<dolfin::Vertex>& self) { self.operator--(); return self; })
      .def("__next__",[](dolfin::MeshEntityIteratorBase<dolfin::Vertex>& self) {
          self.operator++();
          if (self.end())
            throw py::stop_iteration("");
          return *self;
        });

    m.def("vertices", [](dolfin::Mesh& mesh)
          { return dolfin::MeshEntityIteratorBase<dolfin::Vertex>(mesh); });
    m.def("vertices", [](dolfin::MeshEntity& meshentity)
          { return dolfin::MeshEntityIteratorBase<dolfin::Vertex>(meshentity); });

    //--------------------------------------------------------------------------
    // dolfin::MeshFunction class
    py::class_<dolfin::MeshFunction<bool>,
               std::shared_ptr<dolfin::MeshFunction<bool>>>
      (m, "MeshFunction_bool", "DOLFIN MeshFunction object")
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t>())
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t, bool>())
      .def("__getitem__", (const bool& (dolfin::MeshFunction<bool>::*)
                           (std::size_t) const)
           &dolfin::MeshFunction<bool>::operator[])
      .def("__setitem__", [](dolfin::MeshFunction<bool>& self,
                             std::size_t index, bool value)
           { self.operator[](index) = value;})
      .def("__getitem__", (const bool& (dolfin::MeshFunction<bool>::*)
                           (const dolfin::MeshEntity&) const)
           &dolfin::MeshFunction<bool>::operator[])
      .def("__setitem__", [](dolfin::MeshFunction<bool>& self,
                             const dolfin::MeshEntity& index, bool value)
           { self.operator[](index) = value;});

    py::class_<dolfin::MeshFunction<std::size_t>,
               std::shared_ptr<dolfin::MeshFunction<std::size_t>>>
      (m, "MeshFunction_sizet", "DOLFIN MeshFunction object")
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t>())
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t, std::size_t>())
      .def("__getitem__", (const std::size_t& (dolfin::MeshFunction<std::size_t>::*)
                           (std::size_t) const)
           &dolfin::MeshFunction<std::size_t>::operator[])
      .def("__setitem__", [](dolfin::MeshFunction<std::size_t>& self,
                             std::size_t index, std::size_t value)
           { self.operator[](index) = value;})
      .def("__getitem__", (const std::size_t& (dolfin::MeshFunction<std::size_t>::*)
                           (const dolfin::MeshEntity&) const)
           &dolfin::MeshFunction<std::size_t>::operator[])
      .def("__setitem__", [](dolfin::MeshFunction<std::size_t>& self,
                             const dolfin::MeshEntity& index, std::size_t value)
           { self.operator[](index) = value;});

    py::class_<dolfin::MeshFunction<double>,
               std::shared_ptr<dolfin::MeshFunction<double>>>
      (m, "MeshFunction_double", "DOLFIN MeshFunction object")
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t>())
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t, double>())
      .def("__getitem__", (const double& (dolfin::MeshFunction<double>::*)
                           (std::size_t) const)
           &dolfin::MeshFunction<double>::operator[])
      .def("__setitem__", [](dolfin::MeshFunction<double>& self,
                             std::size_t index, double value)
           { self.operator[](index) = value;})
      .def("__getitem__", (const double& (dolfin::MeshFunction<double>::*)
                           (const dolfin::MeshEntity&) const)
           &dolfin::MeshFunction<double>::operator[])
      .def("__setitem__", [](dolfin::MeshFunction<double>& self,
                             const dolfin::MeshEntity& index, double value)
           { self.operator[](index) = value;});

    //--------------------------------------------------------------------------
    // dolfin::MeshEditor class
    py::class_<dolfin::MeshEditor, std::shared_ptr<dolfin::MeshEditor>>
      (m, "MeshEditor", "DOLFIN MeshEditor object")
      .def(py::init<>())
      .def("open", (void (dolfin::MeshEditor::*)(dolfin::Mesh& , std::string, std::size_t, std::size_t, std::size_t))
           &dolfin::MeshEditor::open,
           py::arg("mesh"), py::arg("type"), py::arg("tdim"), py::arg("gdim"), py::arg("degree") = 1)
      .def("init_vertices", &dolfin::MeshEditor::init_vertices)
      .def("init_cells", &dolfin::MeshEditor::init_cells)
      .def("add_vertex", (void (dolfin::MeshEditor::*)(std::size_t, const dolfin::Point&))
           &dolfin::MeshEditor::add_vertex)
      .def("add_cell", (void (dolfin::MeshEditor::*)(std::size_t, const std::vector<std::size_t>&))
           &dolfin::MeshEditor::add_cell)
      .def("close", &dolfin::MeshEditor::close, py::arg("order") = true);


    //--------------------------------------------------------------------------
    // dolfin::SubDomain class

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
           &dolfin::SubDomain::inside);

  }

}
