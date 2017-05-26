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
      _counter(0), _xml_doc(new pugi::xml_document)
{
  add_dolfin_parameters();
  restart_xdmf();
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
void XDMFFile2::restart_xdmf()
{
  // Reset current XML state
  _xml_doc->reset();

  // Create doctype for the XDMF document
  _xml_doc->append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");

  // Add main parent XDMF node
  _xml_doc->append_child("Xdmf");

  // Add version and XInclude
  xdmf_node().append_attribute("Version") = "3.0";
  xdmf_node().append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";

  // Add Domain node
  xdmf_node().append_child("Domain");

}
//-----------------------------------------------------------------------------
pugi::xml_node XDMFFile2::xdmf_node()
{
  return _xml_doc->child("Xdmf");
}
//-----------------------------------------------------------------------------
pugi::xml_node XDMFFile2::domain_node()
{
  return xdmf_node().child("Domain");
}
//----------------------------------------------------------------------------
// FIXME: This should be in HDF5 classes...
std::string XDMFFile2::get_hdf5_filename() const
{
  boost::filesystem::path p(_filename);
  p.replace_extension(".h5");
  if (p.string() == _filename) {
    dolfin_error("XDMFFile.cpp",
                 "deduce name of HDF5 file from XDMF filename",
                 "Filename clash. Check XDMF filename.");
  }

  return p.string();
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

  // Restart the XML file
  restart_xdmf();

  // Retrieve the Domain node
  pugi::xml_node _domain_node = domain_node();

  // Add mesh XML trees to domain node
  pugi::xml_node grid_node = append_grid_node(_domain_node, mesh);
  pugi::xml_node topology_node = append_topology_node(grid_node, mesh);
  pugi::xml_node geometry_node = append_geometry_node(grid_node, mesh);

//  ///
//  /// preparing H5 path and topology data vector
//  ///
//
//  // Compute packed topology data
//  std::vector<T> topology_data;
//
//  // FIXME: bad mesh degree handling
//  dolfin_assert(degree == 1 or degree == 2);
//  if (degree == 1)
//    topology_data = compute_topology_data<T>(mesh, cell_dim);
//  else
//    topology_data = compute_quadratic_topology<T>(mesh);
//
//  // Add topology DataItem node
//  const std::string group_name = path_prefix + "/" + mesh.name();
//  const std::string h5_path = group_name + "/topology";
//  const std::vector<std::int64_t> shape = {num_cells, num_nodes_per_cell};
//  const std::string number_type = "UInt";
//
//  add_data_item(comm, topology_node, h5_id, h5_path,
//                topology_data, shape, number_type);
//
//  ///
//  /// preparing H5 path and geometry data vector
//  ///
//
//  // Pack geometry data
//  std::vector<double> x;
//  if (degree == 1)
//    x = DistributedMeshTools::reorder_vertices_by_global_indices(mesh);
//  else
//    x = mesh_geometry.x();
//
//  // XDMF does not support 1D, so handle as special case
//  if (gdim == 1) {
//    // Pad the coordinates with zeros for a dummy Y
//    gdim = 2;
//    std::vector<double> _x(2 * x.size(), 0.0);
//    for (std::size_t i = 0; i < x.size(); ++i)
//      _x[2 * i] = x[i];
//    std::swap(x, _x);
//  }
//
//  // Add geometry DataItem node
//  const std::string group_name = path_prefix + "/" + mesh.name();
//  const std::string h5_path = group_name + "/geometry";
//  const std::vector<std::int64_t> shape = {num_points, gdim};
//
//  add_data_item(comm, geometry_node, h5_id, h5_path, x, shape);
//
//  ////////

  // Save XML file (on process 0 only)
  if (MPI::rank(_mpi_comm) == 0)
    _xml_doc->save_file(_filename.c_str(), "  ");

}
//----------------------------------------------------------------------------
pugi::xml_node XDMFFile2::append_grid_node(pugi::xml_node &xml_node, const Mesh &mesh)
{
  // Add grid node
  pugi::xml_node grid_node = xml_node.append_child("Grid");

  // Add mesh name and type
  grid_node.append_attribute("Name") = mesh.name().c_str();
  grid_node.append_attribute("GridType") = "Uniform";

  return grid_node;
}
//----------------------------------------------------------------------------
pugi::xml_node XDMFFile2::append_topology_node(pugi::xml_node &xml_node, const Mesh &mesh)
{
  // Get number of cells (global) and vertices per cell from mesh
  // We have to know this to set correct Topology node attributes
  const size_t tdim = mesh.topology().dim();
  const std::int64_t num_cells = mesh.topology().size_global(tdim);
  size_t num_nodes_per_cell = mesh.type().num_vertices(tdim);
  const size_t degree = mesh.geometry().degree();

  // FIXME: this should not be here... this is a logic of Mesh
  if (degree == 2) {
    dolfin_assert(cell_dim == (int) mesh.topology().dim());
    num_nodes_per_cell += mesh.type().num_entities(1);
  }

  // Get VTK string for cell type
  const std::string vtk_cell_str = vtk_cell_type_str(mesh.type().entity_type(tdim), degree);

  // Append with the Topology node
  pugi::xml_node topology_node = xml_node.append_child("Topology");
  dolfin_assert(topology_node);
  topology_node.append_attribute("NumberOfElements") = std::to_string(num_cells).c_str();
  topology_node.append_attribute("TopologyType") = vtk_cell_str.c_str();
  topology_node.append_attribute("NodesPerElement") = num_nodes_per_cell;

  return topology_node;
}
//-----------------------------------------------------------------------------
pugi::xml_node XDMFFile2::append_geometry_node(pugi::xml_node &xml_node, const Mesh &mesh)
{
  const MeshGeometry &mesh_geometry = mesh.geometry();
  size_t gdim = mesh_geometry.dim();

  // Compute number of points (global) in mesh (equal to number of vertices
  // for affine meshes)
  const size_t degree = mesh_geometry.degree();
  dolfin_assert(degree == 1 or degree == 2);
  const std::int64_t num_points
      = (degree == 1) ? mesh.size_global(0) : (mesh.size(0) + mesh.size(1));

  // Add geometry node and attributes
  pugi::xml_node geometry_node = xml_node.append_child("Geometry");
  dolfin_assert(geometry_node);
  dolfin_assert(gdim > 0 and gdim <= 3);
  const std::string geometry_type = (gdim == 3) ? "XYZ" : "XY";
  geometry_node.append_attribute("GeometryType") = geometry_type.c_str();

  return geometry_node;
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
//template<typename T>
//void XDMFFile::append_data_item(pugi::xml_node &xml_node,
//                                hid_t h5_id, const std::string h5_path, const T &x,
//                                const std::vector<std::int64_t> shape,
//                                const std::string number_type)
//{
//  // Add DataItem node
//  dolfin_assert(xml_node);
//  pugi::xml_node data_item_node = xml_node.append_child("DataItem");
//  dolfin_assert(data_item_node);
//
//  // Add dimensions attribute
//  data_item_node.append_attribute("Dimensions")
//      = container_to_string(shape, " ", 16).c_str();
//
//  // Set type for topology data (needed by XDMF to prevent default to float)
//  if (!number_type.empty())
//    data_item_node.append_attribute("NumberType") = number_type.c_str();
//
//  // Add format attribute
//  if (h5_id < 0) {
//    data_item_node.append_attribute("Format") = "XML";
//    dolfin_assert(shape.size() == 2);
//    data_item_node.append_child(pugi::node_pcdata)
//        .set_value(container_to_string(x, " ", 16, shape[1]).c_str());
//  } else {
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
//  }
//}
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