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

#include "PugiXDMFXMLDocument.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
PugiXDMFXMLDocument::PugiXDMFXMLDocument(MPI_Comm comm) : _xml_doc(new pugi::xml_document)
{
}
//-----------------------------------------------------------------------------
void PugiXDMFXMLDocument::reset_xml()
{
  _xml_doc->reset();
}
//-----------------------------------------------------------------------------
void PugiXDMFXMLDocument::add_xdmf_doctype()
{
  _xml_doc->append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");
}
//-----------------------------------------------------------------------------
void PugiXDMFXMLDocument::add_xdmf_node()
{
  pugi::xml_node xdmf_node = _xml_doc->append_child("Xdmf");

  xdmf_node.append_attribute("Version") = "3.0";
  xdmf_node.append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";

}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::get_xdmf_node()
{
  return _xml_doc->child("Xdmf");
}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::add_domain_to_node(pugi::xml_node node, std::string name)
{
  pugi::xml_node domain_node = node.append_child("Domain");
  domain_node.append_attribute("Name") = name.c_str();

  return domain_node;
}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::add_grid_to_node(pugi::xml_node node,
                                                     std::string name,
                                                     std::string grid_type,
                                                     std::string collection_type,
                                                     std::string section)
{
  pugi::xml_node grid_node = node.append_child("Grid");
  grid_node.append_attribute("Name") = name.c_str();

  std::string task = "append grid node to XML document";

  if (_grid_types.find(grid_type) == _grid_types.end()) {
    dolfin_error(__FILE__,
                 task,
                 "GridType \"%s\"",
                 grid_type,
                 "is not supported by XDMF3 standard.");
  }

  if (_collection_types.find(collection_type) == _collection_types.end()) {
    dolfin_error(__FILE__,
                 task,
                 "CollectionType \"%s\"",
                 collection_type,
                 "is not supported by XDMF3 standard.");
  }

  if (_sections.find(section) == _sections.end()) {
    dolfin_error(__FILE__,
                 task,
                 "Section \"%s\"",
                 section,
                 "is not supported by XDMF3 standard.");
  }

  grid_node.append_attribute("GridType") = grid_type.c_str();

  if(grid_type == "Collection")
    grid_node.append_attribute("CollectionType") = collection_type.c_str();

  if(grid_type == "Subset")
    grid_node.append_attribute("Section") = section.c_str();

  return grid_node;
}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::add_topology_to_node(pugi::xml_node node,
                                                         std::string name,
                                                         std::string topology_type,
                                                         std::int8_t number_of_elements,
                                                         std::int8_t nodes_per_element,
                                                         std::string dimensions,
                                                         std::string order)
{
  pugi::xml_node topology_node = node.append_child("Topology");
  topology_node.append_attribute("Name") = name.c_str();

  std::string task = "append topology node to XML document";

  if (_topology_types.find(topology_type) == _topology_types.end()) {
    dolfin_error(__FILE__,
                 task,
                 "TopologyType \"%s\"",
                 topology_type,
                 "is not supported by XDMF3 standard.");
  }

  topology_node.append_attribute("TopologyType") = topology_type.c_str();
  topology_node.append_attribute("NumberOfElements") = std::to_string(number_of_elements).c_str();

  if(!nodes_per_element.empty())
    topology_node.append_attribute("NodesPerElement") = std::to_string(nodes_per_element).c_str();
  if(!dimensions.empty())
    topology_node.append_attribute("Dimensions") = dimensions.c_str();
  if(!order.empty())
    topology_node.append_attribute("Order") = order.c_str();
}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::add_geometry_to_node(pugi::xml_node node,
                                                         std::string geometry_type)
{
  pugi::xml_node geom_node = node.append_child("Geometry");

  std::string task = "append geometry node to XML document";

  if (_geometry_types.find(geometry_type) == _geometry_types.end()) {
    dolfin_error(__FILE__,
                 task,
                 "GeometryType \"%s\"",
                 geometry_type,
                 "is not supported by XDMF3 standard.");
  }

  geom_node.append_attribute("GeometryType") = geometry_type.c_str();
}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::add_data_item_to_node(pugi::xml_node node, std::string name,
                                                          std::string item_type,
                                                          std::string dimensions,
                                                          std::string number_type,
                                                          std::string format,
                                                          std::string precision,
                                                          std::string endian,
                                                          std::string compression,
                                                          std::string seek)
{
  pugi::xml_node data_item_node = node.append_child("DataItem");

  std::string task = "append data item to XML document";

  if (_item_types.find(item_type) == _item_types.end()) {
    dolfin_error(__FILE__,
                 task,
                 "ItemType \"%s\"",
                 item_type,
                 "is not supported by XDMF3 standard.");
  }

  if (_number_types.find(number_type) == _number_types.end()) {
    dolfin_error(__FILE__,
                 task,
                 "NumberType \"%s\"",
                 number_type,
                 "is not supported by XDMF3 standard.");
  }

  if (_formats.find(format) == _formats.end()) {
    dolfin_error(__FILE__,
                 task,
                 "Format \"%s\"",
                 format,
                 "is not supported by XDMF3 standard.");
  }

  if (_endians.find(endian) == _endians.end()) {
    dolfin_error(__FILE__,
                 task,
                 "Endian \"%s\"",
                 endian,
                 "is not supported by XDMF3 standard.");
  }

  if (_compressions.find(compression) == _compressions.end()) {
    dolfin_error(__FILE__,
                 task,
                 "Compression \"%s\"",
                 compression,
                 "is not supported by XDMF3 standard.");
  }

  data_item_node.append_attribute("ItemType") = item_type.c_str();
  data_item_node.append_attribute("Dimensions") = dimensions.c_str();
  data_item_node.append_attribute("NumberType") = number_type.c_str();
  data_item_node.append_attribute("Format") = format.c_str();

  data_item_node.append_attribute("Precision") = precision.c_str();

  // Endian and Compression valid only for Binary Format
  if(format == "Binary") {
    data_item_node.append_attribute("Endian") = endian.c_str();
    data_item_node.append_attribute("Compression") = compression.c_str();
  }

  // Seek valid only for Binary Format and Raw Compression
  if(format == "Binary" and compression == "Raw")
    data_item_node.append_attribute("Seek") = seek.c_str();

  return data_item_node;
}
//-----------------------------------------------------------------------------
pugi::xml_node PugiXDMFXMLDocument::add_time_to_node(pugi::xml_node node,
                                                     std::string time_type,
                                                     std::string value)
{
  pugi::xml_node time_node = node.append_child("Time");

  if (_time_types.find(time_type) == _time_types.end()) {
    dolfin_error(__FILE__,
                 "append time item node to XML document",
                 "TimeType \"%s\"",
                 time_type,
                 "is not supported by XDMF3 standard.");
  }

  time_node.append_attribute("TimeType") = time_type.c_str();

  // Value valid only for Single TimeType
  if(time_type == "Single")
    time_node.append_attribute("Value") = value.c_str();

  return time_node;
}
//-----------------------------------------------------------------------------
void PugiXDMFXMLParser::save_file(std::string xml_filename)
{
  _xml_doc->save_file(xml_filename.c_str(), "  ");
}
//-----------------------------------------------------------------------------