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

#include <dolfin/common/Variable.h>
#include <dolfin/geometry/BoundingBoxTree.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MultiMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include <dolfin/mesh/MeshQuality.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/SubMesh.h>

#include "../mpi_interface.h"

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
    //-------------------------------------------------------------------------

    // dolfin::Mesh class
    py::class_<dolfin::Mesh, std::shared_ptr<dolfin::Mesh>>(m, "Mesh", py::dynamic_attr(), "DOLFIN Mesh object")
      .def(py::init<>())
      .def("bounding_box_tree", &dolfin::Mesh::bounding_box_tree)
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
             return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
               (self.geometry().x().data(),
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
      .def("num_edges", &dolfin::Mesh::num_edges, "Number of edges")
      .def("num_faces", &dolfin::Mesh::num_faces, "Number of faces")
      .def("num_facets", &dolfin::Mesh::num_facets, "Number of facets")
      .def("num_cells", &dolfin::Mesh::num_cells, "Number of cells")
      .def("rmax", &dolfin::Mesh::rmax)
      .def("rmin", &dolfin::Mesh::rmin)
      .def("size_global", &dolfin::Mesh::size_global)
      .def("topology", (const dolfin::MeshTopology& (dolfin::Mesh::*)() const)
           &dolfin::Mesh::topology, "Mesh topology")
      .def("translate", &dolfin::Mesh::translate)
      // UFL related
      .def("ufl_id", [](const dolfin::Mesh& self){ return self.id(); })
      .def("cell_name", [](const dolfin::Mesh& self)
           { return dolfin::CellType::type2string(self.type().cell_type()); }
        );

    //-------------------------------------------------------------------------
    // dolfin::BoundaryMesh class
    py::class_<dolfin::BoundaryMesh, std::shared_ptr<dolfin::BoundaryMesh>, dolfin::Mesh>
      (m, "BoundaryMesh", "DOLFIN BoundaryMesh object")
      .def(py::init<const dolfin::Mesh&, std::string, bool>(),
           py::arg("mesh"), py::arg("type"), py::arg("order")=true);

    //-------------------------------------------------------------------------
    // dolfin::MeshConnectivity class
    py::class_<dolfin::MeshConnectivity, std::shared_ptr<dolfin::MeshConnectivity>>
      (m, "MeshConnectivity", "DOLFIN MeshConnectivity object");

    //-------------------------------------------------------------------------
    // dolfin::MeshTopology class
    py::class_<dolfin::MeshTopology, std::shared_ptr<dolfin::MeshTopology>>
      (m, "MeshTopology", "DOLFIN MeshTopology object")
      .def("dim", &dolfin::MeshTopology::dim, "Topological dimension")
      .def("__call__", (const dolfin::MeshConnectivity& (dolfin::MeshTopology::*)(std::size_t, std::size_t) const)
           &dolfin::MeshTopology::operator());

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
      .def("dim", &dolfin::MeshEntity::dim, "Topological dimension")
      .def("mesh", &dolfin::MeshEntity::mesh, "Associated mesh")
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
      .def(py::init<const dolfin::Mesh&, std::size_t>())
      .def("point", &dolfin::Vertex::point);

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
      .def("collides", (bool (dolfin::Cell::*)(const dolfin::Point&) const) &dolfin::Cell::collides)
      .def("collides", (bool (dolfin::Cell::*)(const dolfin::MeshEntity&) const) &dolfin::Cell::collides)
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

    m.def("entities", [](dolfin::Mesh& mesh, std::size_t dim)
          { return dolfin::MeshEntityIterator(mesh, dim); });
    m.def("entities", [](dolfin::MeshEntity& meshentity, std::size_t dim)
          { return dolfin::MeshEntityIterator(meshentity, dim); });

#define MESHITERATOR_MACRO(TYPE, NAME) \
    py::class_<dolfin::MeshEntityIteratorBase<dolfin::TYPE>, \
               std::shared_ptr<dolfin::MeshEntityIteratorBase<dolfin::TYPE>>> \
      (m, #TYPE"Iterator", "DOLFIN "#TYPE"Iterator object") \
      .def(py::init<const dolfin::Mesh&>()) \
      .def("__iter__",[](dolfin::MeshEntityIteratorBase<dolfin::TYPE>& self) { self.operator--(); return self; }) \
      .def("__next__",[](dolfin::MeshEntityIteratorBase<dolfin::TYPE>& self) { \
          self.operator++(); \
          if (self.end()) \
            throw py::stop_iteration(""); \
          return *self; \
        }); \
 \
    m.def(#NAME, [](dolfin::Mesh& mesh) \
          { return dolfin::MeshEntityIteratorBase<dolfin::TYPE>(mesh); }); \
    m.def(#NAME, [](dolfin::MeshEntity& meshentity)\
          { return dolfin::MeshEntityIteratorBase<dolfin::TYPE>(meshentity); })

    MESHITERATOR_MACRO(Cell, cells);
    MESHITERATOR_MACRO(Facet, facets);
    MESHITERATOR_MACRO(Face, faces);
    MESHITERATOR_MACRO(Edge, edges);
    MESHITERATOR_MACRO(Vertex, vertices);
#undef MESHITERATOR_MACRO

    //--------------------------------------------------------------------------
    // dolfin::MeshFunction class

#define MESHFUNCTION_MACRO(SCALAR, SCALAR_NAME) \
    py::class_<dolfin::MeshFunction<SCALAR>, \
      std::shared_ptr<dolfin::MeshFunction<SCALAR>>, dolfin::Variable>  \
      (m, "MeshFunction_"#SCALAR_NAME, "DOLFIN MeshFunction object") \
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t>()) \
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t, SCALAR>()) \
      .def("__getitem__", (const SCALAR& (dolfin::MeshFunction<SCALAR>::*) \
                           (std::size_t) const) \
           &dolfin::MeshFunction<SCALAR>::operator[]) \
      .def("__setitem__", [](dolfin::MeshFunction<SCALAR>& self, \
                             std::size_t index, SCALAR value) \
           { self.operator[](index) = value;}) \
      .def("__getitem__", (const SCALAR& (dolfin::MeshFunction<SCALAR>::*) \
                           (const dolfin::MeshEntity&) const) \
           &dolfin::MeshFunction<SCALAR>::operator[]) \
      .def("__setitem__", [](dolfin::MeshFunction<SCALAR>& self, \
                             const dolfin::MeshEntity& index, SCALAR value) \
           { self.operator[](index) = value;}) \
      .def("__len__", &dolfin::MeshFunction<SCALAR>::size) \
      .def("size", &dolfin::MeshFunction<SCALAR>::size) \
      .def("set_all", &dolfin::MeshFunction<SCALAR>::set_all) \
      .def("where_equal", &dolfin::MeshFunction<SCALAR>::where_equal) \
      .def("array", [](dolfin::MeshFunction<SCALAR>& self) \
           { return Eigen::Map<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>(self.values(), self.size()); })

    MESHFUNCTION_MACRO(bool, bool);
    MESHFUNCTION_MACRO(int, int);
    MESHFUNCTION_MACRO(double, double);
    MESHFUNCTION_MACRO(std::size_t, sizet);
#undef MESHFUNCTION_MACRO

#define MESH_ENTITY_FUNCTION_MACRO(TYPE, SCALAR, SCALAR_NAME) \
    py::class_<dolfin::TYPE<SCALAR>, std::shared_ptr<dolfin::TYPE<SCALAR>>, \
      dolfin::MeshFunction<SCALAR>>(m, #TYPE"_"#SCALAR_NAME)

    MESH_ENTITY_FUNCTION_MACRO(VertexFunction, bool, bool);
    MESH_ENTITY_FUNCTION_MACRO(VertexFunction, int, int);
    MESH_ENTITY_FUNCTION_MACRO(VertexFunction, double, double);
    MESH_ENTITY_FUNCTION_MACRO(VertexFunction, std::size_t, sizet);
    MESH_ENTITY_FUNCTION_MACRO(EdgeFunction, bool, bool);
    MESH_ENTITY_FUNCTION_MACRO(EdgeFunction, int, int);
    MESH_ENTITY_FUNCTION_MACRO(EdgeFunction, double, double);
    MESH_ENTITY_FUNCTION_MACRO(EdgeFunction, std::size_t, sizet);
    MESH_ENTITY_FUNCTION_MACRO(FaceFunction, bool, bool);
    MESH_ENTITY_FUNCTION_MACRO(FaceFunction, int, int);
    MESH_ENTITY_FUNCTION_MACRO(FaceFunction, double, double);
    MESH_ENTITY_FUNCTION_MACRO(FaceFunction, std::size_t, sizet);
    MESH_ENTITY_FUNCTION_MACRO(FacetFunction, bool, bool);
    MESH_ENTITY_FUNCTION_MACRO(FacetFunction, int, int);
    MESH_ENTITY_FUNCTION_MACRO(FacetFunction, double, double);
    MESH_ENTITY_FUNCTION_MACRO(FacetFunction, std::size_t, sizet);
    MESH_ENTITY_FUNCTION_MACRO(CellFunction, bool, bool);
    MESH_ENTITY_FUNCTION_MACRO(CellFunction, int, int);
    MESH_ENTITY_FUNCTION_MACRO(CellFunction, double, double);
    MESH_ENTITY_FUNCTION_MACRO(CellFunction, std::size_t, sizet);
#undef MESH_ENTITY_FUNCTION_MACRO

    //--------------------------------------------------------------------------
    // dolfin::MeshValueCollection class
#define MESHVALUECOLLECTION_MACRO(SCALAR, SCALAR_NAME) \
    py::class_<dolfin::MeshValueCollection<SCALAR>, \
      std::shared_ptr<dolfin::MeshValueCollection<SCALAR>>, dolfin::Variable>   \
      (m, "MeshValueCollection_"#SCALAR_NAME, "DOLFIN MeshValueCollection object") \
      .def(py::init<std::shared_ptr<const dolfin::Mesh>>()) \
      .def(py::init<std::shared_ptr<const dolfin::Mesh>, std::size_t>()) \
      .def("dim", &dolfin::MeshValueCollection<SCALAR>::dim) \
      .def("size", &dolfin::MeshValueCollection<SCALAR>::size) \
      .def("get_value", &dolfin::MeshValueCollection<SCALAR>::get_value) \
      .def("set_value", (bool (dolfin::MeshValueCollection<SCALAR>::*)(std::size_t, const SCALAR&)) \
           &dolfin::MeshValueCollection<SCALAR>::set_value) \
      .def("assign", [](dolfin::MeshValueCollection<SCALAR>& self, const dolfin::MeshValueCollection<SCALAR>& other) \
         { self = other; })

    MESHVALUECOLLECTION_MACRO(bool, bool);
    MESHVALUECOLLECTION_MACRO(int, int);
    MESHVALUECOLLECTION_MACRO(double, double);
    MESHVALUECOLLECTION_MACRO(std::size_t, sizet);
#undef MESHVALUECOLLECTION_MACRO

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
      .def("add_vertex", (void (dolfin::MeshEditor::*)(std::size_t, const std::vector<double>&))
           &dolfin::MeshEditor::add_vertex)
      .def("add_cell", (void (dolfin::MeshEditor::*)(std::size_t, const std::vector<std::size_t>&))
           &dolfin::MeshEditor::add_cell)
      .def("close", &dolfin::MeshEditor::close, py::arg("order") = true);

    //--------------------------------------------------------------------------
    // dolfin::MultiMesh class
    py::class_<dolfin::MultiMesh, std::shared_ptr<dolfin::MultiMesh>>
      (m, "MultiMesh", "DOLFIN MultiMesh")
      .def(py::init<>());

    //--------------------------------------------------------------------------
    // dolfin::MeshQuality class
    py::class_<dolfin::MeshQuality>
      (m, "MeshQuality", "DOLFIN MeshQuality class")
      .def_static("radius_ratios", &dolfin::MeshQuality::radius_ratios)
      .def_static("radius_ratio_histogram_data", &dolfin::MeshQuality::radius_ratio_histogram_data)
      .def_static("radius_ratio_min_max", &dolfin::MeshQuality::radius_ratio_min_max)
      .def_static("radius_ratio_matplotlib_histogram", &dolfin::MeshQuality::radius_ratio_matplotlib_histogram)
      .def_static("dihedral_angles_min_max", &dolfin::MeshQuality::dihedral_angles_min_max)
      .def_static("dihedral_angles_matplotlib_histogram", &dolfin::MeshQuality::dihedral_angles_matplotlib_histogram);

    //--------------------------------------------------------------------------
    // dolfin::SubMesh class
    py::class_<dolfin::SubMesh, std::shared_ptr<dolfin::SubMesh>>
      (m, "SubMesh", "DOLFIN SubMesh")
      .def(py::init<const dolfin::Mesh&, const dolfin::SubDomain&>());

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
