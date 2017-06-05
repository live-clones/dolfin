// Copyright (C) 2017 M. Habera
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
#ifndef DOLFIN_PUGIXDMFXMLDOCUMENT_H
#define DOLFIN_PUGIXDMFXMLDOCUMENT_H

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dolfin/common/MPI.h>
#include "pugixml.hpp"

namespace dolfin
{

/// Prepares and saves XDMF XML skeleton document
///
/// This class wraps the functionality of Pugi XML library
/// http://pugixml.org/ and defines standards for XDMF3 document.
///
/// Current standards for XDMF3 document are taken from
/// http://www.xdmf.org/index.php/XDMF_Model_and_Format
/// and are stored as private member variables.
///
/// This class is responsible for building valid XML XDMF document.
///

class PugiXDMFXMLDocument
{
public:

  PugiXDMFXMLDocument(MPI_Comm comm);

  /// Resets underlying Pugi XML document
  void reset_xml();

  /// Adds XDMF doctype tag to Pugi XML document
  ///
  /// Appends root of XML document with <!DOCTYPE> node.
  void add_xdmf_doctype();

  /// Adds XDMF tag to Pugi XML document
  ///
  /// Appends root of XML document with <Xdmf/> node
  /// and sets the xdmf version=3.0 and XInclude namespace from
  /// [http://www.w3.org/2001/XInclude]
  void add_xdmf_node();

  /// Returns Xdmf node
  pugi::xml_node get_xdmf_node();

  /// Adds Domain node to a given node
  ///
  /// There could possibly be several domains per each XDMF file.
  ///
  /// @param node A node to which Domain node is appended
  /// @param name Name attribute of Domain node
  /// @return Newly appended Domain node
  pugi::xml_node add_domain_to_node(pugi::xml_node node,
                                    std::string name);

  /// Adds Grid node to a given node
  ///
  /// A Grid is a container for information related to 2D and 3D points,
  /// structured or unstructured connectivity, and assigned values.
  ///
  /// @param node A node to which Grid node is appended
  /// @param name Name attribute of Grid node
  /// @param grid_type GridType attribute value
  /// @param collection_type CollectionType attribute value. Ignored if GridType != Collection.
  /// @param section Section attribute value. Ignored if GridType != SubSet.
  /// @return Newly appended Grid node
  pugi::xml_node add_grid_to_node(pugi::xml_node node,
                                  std::string name,
                                  std::string grid_type,
                                  std::string collection_type = "",
                                  std::string section = "");

  /// Adds Topology node to a given node
  ///
  /// @param node A node to which Topology node is appended
  /// @param name Name attribute of Topology node
  /// @param topology_type TopologyType attribute value
  /// @param nodes_per_element NodesPerElement attribute value
  /// @param number_of_elements NumberOfElements attribute value
  /// @param dimensions Dimensions attribute value
  /// @param order Order attribute value
  /// @return Newly appended Topology node
  pugi::xml_node add_topology_to_node(pugi::xml_node node,
                                      std::string name,
                                      std::string topology_type,
                                      std::int64_t number_of_elements,
                                      std::int8_t nodes_per_element = 0,
                                      std::string dimensions = "",
                                      std::string order = "");

  /// Adds Geometry node to a given node
  ///
  /// @param node A node to which Geometry node is appended
  /// @param geometry_type GeometryType attribute value
  /// @return geometry_node Newly appended geometry node
  pugi::xml_node add_geometry_to_node(pugi::xml_node node,
                                      std::string geometry_type);

  /// Adds DataItem node to a given node
  ///
  /// @param node A node to which DataItem node is appended
  /// @param name Name attribute value
  /// @param item_type ItemType attribute value
  /// @param dimensions Dimensions attribute value
  /// @param number_type NumberType attribute value
  /// @param format Format attribute value
  /// @param precision Precision attribute value
  /// @param endian Endian attribute value. Ignored if Format != Binary.
  /// @param compression Compression attribute value. Ignored if Format != Binary.
  /// @param seek Seek attribute value. Ignored if Format != Binary and Compression != Raw.
  /// @return Newly appended DataItem node
  pugi::xml_node add_data_item_to_node(pugi::xml_node node,
                                       std::string name,
                                       std::string item_type,
                                       std::string dimensions,
                                       std::string number_type,
                                       std::string format,
                                       std::string precision,
                                       std::string endian = "",
                                       std::string compression = "",
                                       std::string seek = "");

  /// Adds Time node to a given node
  ///
  /// @param node A node to which Time node is appended
  /// @param time_type TimeType attribute value
  /// @param value Value attribute value. Ignored if TimeType != Single.
  /// @return Newly appended Time node
  pugi::xml_node add_time_to_node(pugi::xml_node node,
                                  std::string time_type,
                                  std::string value);

  void set_data(pugi::xml_node node,
                std::string data);

  /// Saves XDMF XML to a file of given filename.
  /// Uses internal Pugi's xml_doc->save_file() method.
  ///
  /// @param xml_filename Filename of XDMF XML file
  void save_file(std::string xml_filename);

private:

  std::unique_ptr<pugi::xml_document> _xml_doc;

  const std::set<std::string> _grid_types = {"Uniform", "Collection", "Tree", "Subset"};
  const std::set<std::string> _collection_types = {"Spatial", "Temporal", ""};
  const std::set<std::string> _sections = {"DataItem", "All", ""};

  const std::set<std::string> _topology_types =
      {"Polyvertex", "Polyline", "Polygon", "Triangle", "Quadrilateral", "Tetrahedron", "Pyramid", "Wedge",
       "Hexahedron", "Edge_3", "Triangle_6", "Quadrilateral_8", "Tetrahedron_10", "Pyramid_13", "Wedge_15",
       "Hexahedron_20", "Mixed", "2DSMesh", "2DRectMesh", "2DCoRectMesh", "3DSMesh", "3DRectMesh", "3DCoRectMesh"};

  const std::set<std::string> _item_types = {"Uniform", "Collection", "tree", "HyperSlab", "coordinates", "Function"};
  const std::set<std::string> _number_types = {"Float", "Int", "UInt", "Char", "UChar"};
  const std::set<std::string> _formats = {"XML", "HDF", "Binary"};
  const std::set<std::string> _endians = {"Native", "Big", "Little", ""};
  const std::set<std::string> _compressions = {"Raw", "Zlib", "BZip2", ""};

  const std::set<std::string> _geometry_types = {"XYZ", "XY", "X_Y_Z", "VxVyVz", "Origin_DxDyDz", "Origin_DxDy"};

  const std::set<std::string> _time_types = {"Single", "HyperSlab", "List", "Range"};

};

}
#endif //DOLFIN_PUGIXDMFXMLDOCUMENT_H
