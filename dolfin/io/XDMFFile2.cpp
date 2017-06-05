// Copyright (C) 2012-2017 Chris N. Richardson, Garth N. Wells and M. Habera
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
//

#include <iomanip>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/container/vector.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include "pugixml.hpp"

#include <dolfin/common/MPI.h>
#include <dolfin/common/defines.h>
#include <dolfin/common/utils.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/DistributedMeshTools.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/LocalMeshData.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/GlobalParameters.h>
#include "HDF5File.h"
#include "HDF5Utility.h"
#include "XDMFFile2.h"
#include "xmlutils.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
XDMFFile2::XDMFFile2(MPI_Comm comm, const std::string filename)
    : _mpi_comm(comm), _filename(filename),
      _counter(0), _xdmf_xml_doc(new PugiXDMFXMLDocument(comm))
{
  add_dolfin_parameters();
}
//-----------------------------------------------------------------------------
XDMFFile2::~XDMFFile2()
{
  close();
}
//-----------------------------------------------------------------------------
void XDMFFile2::close()
{
  if (has_hdf5()) {
    // Close the HDF5 file
    _hdf5_file.reset();
  }
}
//-----------------------------------------------------------------------------
void XDMFFile2::add_dolfin_parameters()
{
  // Rewrite the mesh at every time step in a time series. Should be
  // turned off if the mesh remains constant.
  parameters.add("rewrite_function_mesh", true);

  // Functions share the same mesh for the same time step. The files
  // produced are smaller and work better in Paraview
  parameters.add("functions_share_mesh", false);

  // FIXME: This is only relevant to HDF5
  // Flush datasets to disk at each timestep. Allows inspection of the
  // HDF5 file whilst running, at some performance cost.
  parameters.add("flush_output", false);
}
//-----------------------------------------------------------------------------
void XDMFFile2::write(const Mesh &mesh, Encoding encoding)
{

  info("Writing mesh to XDMF file in %s ...",
       MPI::size(_mpi_comm) == 1 ? "serial" : "parallel, rank = %i",
       MPI::rank(_mpi_comm));

  // Check that encoding is supported
  check_encoding(encoding);

  // Check that mesh is supported
  check_mesh(mesh);

  ///
  ///  get some info from mesh needed for XML XDMF
  ///

  // Get number of cells (global) and vertices per cell from mesh
  // We have to know this to set correct Topology node attributes
  const size_t tdim = mesh.topology().dim();
  size_t gdim = mesh.geometry().dim();
  const std::int64_t num_cells = mesh.topology().size_global(tdim);
  size_t num_nodes_per_cell = mesh.type().num_vertices(tdim);
  const size_t degree = mesh.geometry().degree();
  // Get VTK string for cell type
  const std::string vtk_cell_str = vtk_cell_type_str(mesh.type().entity_type(tdim), degree);
  const std::string geometry_type = (gdim == 3) ? "XYZ" : "XY";

  // FIXME: this should not be here... this is a logic of Mesh
  if (degree == 2) {
    num_nodes_per_cell += mesh.type().num_entities(1);
  }

  // FIXME: allow arbitrary mesh degree, no this hardcoded stuff...
  // Compute number of points (global) in mesh (equal to number of vertices
  // for affine meshes)
  const std::int64_t num_points = (degree == 1) ? mesh.size_global(0) : (mesh.size(0) + mesh.size(1));

//  const std::string h5_path = "Mesh/" + mesh.name() + "/topology";
  std::vector<std::int64_t> shape = {num_cells, num_nodes_per_cell};
  const std::string dimensions = container_to_string(shape, " ", 16).c_str();

  std::string format;
  if (encoding == Encoding::ASCII)
    format = "XML";
  else if (encoding == Encoding::HDF5)
    format = "HDF";

  ///
  ///  prepare XDMF XML
  ///

  // Reset current XML state
  _xdmf_xml_doc->reset_xml();

  // Create doctype for the XDMF document
  _xdmf_xml_doc->add_xdmf_doctype();

  // Add main parent XDMF node
  _xdmf_xml_doc->add_xdmf_node();

  pugi::xml_node xdmf_node = _xdmf_xml_doc->get_xdmf_node();
  pugi::xml_node domain_node = _xdmf_xml_doc->add_domain_to_node(xdmf_node, "Domain");
  pugi::xml_node grid_node = _xdmf_xml_doc->add_grid_to_node(domain_node, "", "Uniform");
  pugi::xml_node topology_node = _xdmf_xml_doc->add_topology_to_node(grid_node,
                                                                     "",
                                                                     vtk_cell_str.c_str(),
                                                                     num_cells);
  pugi::xml_node geometry_node = _xdmf_xml_doc->add_geometry_to_node(grid_node, geometry_type);

  pugi::xml_node topology_data_item_node = _xdmf_xml_doc->add_data_item_to_node(topology_node, "", "Uniform",
                                                                                dimensions, "UInt", format, "4");
  pugi::xml_node geometry_data_item_node = _xdmf_xml_doc->add_data_item_to_node(geometry_node, "", "Uniform",
                                                                                dimensions, "Float", format, "4");

  ///
  /// preparing topology data vector
  ///

  // Compute packed topology data
  std::vector<std::int64_t> topology_data;

  // FIXME: bad mesh degree handling
  if (degree == 1) {

    topology_data.reserve(mesh.num_entities(tdim) * (num_nodes_per_cell));
    const std::vector<std::int8_t> perm = mesh.type().vtk_mapping();

    if (MPI::size(_mpi_comm) == 1) {
      // Simple case when nothing is shared between processes
      if (tdim == 0) {
        for (VertexIterator v(mesh); !v.end(); ++v)
          topology_data.push_back(v->global_index());
      } else {
        const auto &global_vertices = mesh.topology().global_indices(0);
        for (MeshEntityIterator c(mesh, tdim); !c.end(); ++c) {
          const unsigned int *entities = c->entities(0);
          for (unsigned int i = 0; i != c->num_entities(0); ++i)
            topology_data.push_back(global_vertices[entities[perm[i]]]);
        }
      }
    } else {
      std::set<unsigned int> non_local_entities = compute_nonlocal_entities(mesh, tdim);

      if (tdim == 0) {
        // Special case for mesh of points
        for (VertexIterator v(mesh); !v.end(); ++v) {
          if (non_local_entities.find(v->index()) == non_local_entities.end())
            topology_data.push_back(v->global_index());
        }
      } else {
        // Local-to-global map for point indices
        const auto &global_vertices = mesh.topology().global_indices(0);
        for (MeshEntityIterator e(mesh, tdim); !e.end(); ++e) {
          // If not excluded, add to topology
          if (non_local_entities.find(e->index()) == non_local_entities.end()) {
            for (unsigned int i = 0; i != e->num_entities(0); ++i) {
              const unsigned int local_idx = e->entities(0)[perm[i]];
              topology_data.push_back(global_vertices[local_idx]);
            }
          }

        }

      }
    }
  } else {
    const MeshGeometry &geom = mesh.geometry();

    if (geom.degree() != 2 or MPI::size(mesh.mpi_comm()) != 1) {
      dolfin_error("XDMFFile.cpp",
                   "create topology data",
                   "XDMF quadratic mesh only supported in serial");
    }

    std::vector<std::size_t> edge_mapping;
    if (tdim == 1)
      edge_mapping = {0};
    else if (tdim == 2)
      edge_mapping = {2, 0, 1};
    else
      edge_mapping = {5, 2, 4, 3, 1, 0};

    // Get number of points per cell
    const CellType &celltype = mesh.type();
    std::size_t npoint = celltype.num_entities(0) + celltype.num_entities(1);
    topology_data.reserve(npoint * mesh.size(tdim));

    for (CellIterator c(mesh); !c.end(); ++c) {
      // Add indices for vertices and edges
      for (unsigned int dim = 0; dim != 2; ++dim) {
        for (unsigned int i = 0; i != celltype.num_entities(dim); ++i) {
          std::size_t im = (dim == 0) ? i : edge_mapping[i];
          const std::size_t entity_index
              = (dim == tdim) ? c->index() : c->entities(dim)[im];
          const std::size_t local_idx = geom.get_entity_index(dim, 0, entity_index);
          topology_data.push_back(local_idx);
        }
      }
    }
  }
//
//

  ///
  ///  Writing topology data vector
  ///

  if (encoding == Encoding::ASCII) {
    _xdmf_xml_doc->set_data(topology_data_item_node, container_to_string(topology_data, " ", 16, shape[1]));
  } else {
//#ifdef HAS_HDF5
//    data_item_node.append_attribute("Format") = "HDF";
//
//    // Get name of HDF5 file
//    const std::string hdf5_filename = HDF5Interface::get_filename(h5_id);
//    const boost::filesystem::path p(hdf5_filename);
//
//    // Add HDF5 filename and HDF5 internal path to XML file
//    const std::string xdmf_path = p.filename().string() + ":" + h5_path;
//    data_item_node.append_child(pugi::node_pcdata).set_value(xdmf_path.c_str());
//
//    // Compute total number of items and check for consistency with shape
//    dolfin_assert(!shape.empty());
//    std::int64_t num_items_total = 1;
//    for (auto n : shape)
//      num_items_total *= n;
//
//    dolfin_assert(num_items_total == (std::int64_t) MPI::sum(comm, x.size()));
//
//    // Compute data offset and range of values
//    std::int64_t local_shape0 = x.size();
//    for (std::size_t i = 1; i < shape.size(); ++i)
//    {
//      dolfin_assert(local_shape0 % shape[i] == 0);
//      local_shape0 /= shape[i];
//    }
//    const std::int64_t offset = MPI::global_offset(comm, local_shape0, true);
//    const std::pair<std::int64_t, std::int64_t> local_range
//      = {offset, offset + local_shape0};
//
//    const bool use_mpi_io = (MPI::size(comm) > 1);
//    HDF5Interface::write_dataset(h5_id, h5_path, x, local_range, shape, use_mpi_io,
//                                 false);
//
//    // Add partitioning attribute to dataset
//    std::vector<std::size_t> partitions;
//    std::vector<std::size_t> offset_tmp(1, offset);
//    MPI::gather(comm, offset_tmp, partitions);
//    MPI::broadcast(comm, partitions);
//    HDF5Interface::add_attribute(h5_id, h5_path, "partition", partitions);
//
//#else
//    // Should never reach this point
//    dolfin_error("XDMFFile.cpp",
//                 "add dataitem",
//                 "DOLFIN has not been configured with HDF5");
//#endif
  }
//
//  add_data_item(comm, topology_node, h5_id, h5_path,
//                topology_data, shape, number_type);
//


  ///
  /// preparing geometry data vector
  ///

  // Pack geometry data
  std::vector<double> geometry_data;
  if (degree == 1)
    geometry_data = DistributedMeshTools::reorder_vertices_by_global_indices(mesh);
  else
    geometry_data = mesh.geometry().x();

  // XDMF does not support 1D, so handle as special case
  if (gdim == 1) {
    // Pad the coordinates with zeros for a dummy Y
    gdim = 2;
    std::vector<double> _x(2 * geometry_data.size(), 0.0);
    for (std::size_t i = 0; i < geometry_data.size(); ++i)
      _x[2 * i] = geometry_data[i];
    std::swap(geometry_data, _x);
  }

  ///
  /// writing geom data vector
  ///

  // Add format attribute
  if (encoding == Encoding::ASCII) {
    _xdmf_xml_doc->set_data(geometry_data_item_node, container_to_string(geometry_data,
                                                                         " ",
                                                                         16,
                                                                         shape[1]).c_str());
  }

  ///
  /// writing XML file
  ///

  _xdmf_xml_doc->save_file(_filename.c_str());

}
//----------------------------------------------------------------------------
std::set<unsigned int> XDMFFile2::compute_nonlocal_entities(const Mesh &mesh, int cell_dim)
{
  // If not already numbered, number entities of
  // order cell_dim so we can get shared_entities
  DistributedMeshTools::number_entities(mesh, cell_dim);

  const std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());
  const std::map<std::int32_t, std::set<unsigned int>> &shared_entities
      = mesh.topology().shared_entities(cell_dim);

  std::set<unsigned int> non_local_entities;

  const std::size_t tdim = mesh.topology().dim();
  bool ghosted
      = (mesh.topology().size(tdim) > mesh.topology().ghost_offset(tdim));

  if (!ghosted) {
    // No ghost cells - exclude shared entities
    // which are on lower rank processes
    for (const auto &e : shared_entities) {
      const unsigned int lowest_rank_owner = *(e.second.begin());
      if (lowest_rank_owner < mpi_rank)
        non_local_entities.insert(e.first);
    }
  } else {
    // Iterate through ghost cells, adding non-ghost entities
    // which are in lower rank process cells
    for (MeshEntityIterator c(mesh, tdim, "ghost"); !c.end(); ++c) {
      const unsigned int cell_owner = c->owner();
      for (MeshEntityIterator e(*c, cell_dim); !e.end(); ++e)
        if (!e->is_ghost() && cell_owner < mpi_rank)
          non_local_entities.insert(e->index());
    }
  }
  return non_local_entities;
}
//----------------------------------------------------------------------------
void XDMFFile2::check_encoding(Encoding encoding) const
{
  if (encoding == Encoding::HDF5 and !has_hdf5()) {
    dolfin_error("XDMFFile.cpp",
                 "write/read XDMF file",
                 "DOLFIN has not been compiled with HDF5 support");
  }

  if (encoding == Encoding::ASCII and MPI::size(_mpi_comm) != 1) {
    dolfin_error("XDMFFile.cpp",
                 "write/read XDMF file",
                 "ASCII format is not supported in parallel, use HDF5");
  }
}
//----------------------------------------------------------------------------
void XDMFFile2::check_mesh(const Mesh &mesh) const
{
  if (mesh.num_cells() == 0) {
    dolfin_error("XDMFFile2.cpp",
                 "write/read mesh to/from XDMF file",
                 "Supplied mesh is empty");
  }

  if (mesh.geometry().degree() > 2) {
    dolfin_error("XDMFFile2.cpp",
                 "write/read mesh to/from XDMF file",
                 "Supplied mesh has degree %i which is not supported in XDMF IO");
  }
}
//-----------------------------------------------------------------------------
std::string XDMFFile2::vtk_cell_type_str(CellType::Type cell_type, const size_t degree)
{
  // FIXME: Move to CellType?
  switch (cell_type) {
    case CellType::Type::point:
      switch (degree) {
        case 1:return "PolyVertex";
      }
    case CellType::Type::interval:
      switch (degree) {
        case 1:return "PolyLine";
        case 2:return "Edge_3";
      }
    case CellType::Type::triangle:
      switch (degree) {
        case 1:return "Triangle";
        case 2:return "Tri_6";
      }
    case CellType::Type::quadrilateral:
      switch (degree) {
        case 1:return "Quadrilateral";
        case 2:return "Quad_8";
      }
    case CellType::Type::tetrahedron:
      switch (degree) {
        case 1:return "Tetrahedron";
        case 2:return "Tet_10";
      }
    case CellType::Type::hexahedron:
      switch (degree) {
        case 1:return "Hexahedron";
        case 2:return "Hex_20";
      }
    default:
      dolfin_error("XDMFFile.cpp",
                   "output mesh topology",
                   "Invalid combination of cell type and order");
      return "error";
  }
}