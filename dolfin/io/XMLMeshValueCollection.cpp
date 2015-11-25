// Copyright (C) 2015 Chris Richardson
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

#include "XMLMeshValueCollection.h"

namespace dolfin {


  template <>
  void
  XMLMeshValueCollection::read(MeshValueCollection<std::complex<double>>& mesh_value_collection,
                               const std::string type,
                               const pugi::xml_node xml_node)
  {
    // Get node
    const pugi::xml_node mvc_node
      = get_node(xml_node, "mesh_value_collection");
    dolfin_assert(mvc_node);

    // Get attributes
    const std::string name = mvc_node.attribute("name").value();
    const std::string type_file = mvc_node.attribute("type").value();
    const std::size_t dim = mvc_node.attribute("dim").as_uint();

    // Attach name to mesh value collection object
    mesh_value_collection.rename(name, "a mesh value collection");

    // Set dimension
    mesh_value_collection.init(dim);

    // Check that types match
    if (type != type_file)
    {
      dolfin_error("XMLMeshValueCollection.h",
                   "read mesh value collection from XML file",
                   "Type mismatch, found \"%s\" but expecting \"%s\"",
                   type_file.c_str(), type.c_str());
    }

    // Clear old values
    mesh_value_collection.clear();

    if (type == "complex")
    {
      pugi::xml_node_iterator it;
      for (it = mvc_node.begin(); it != mvc_node.end(); ++it)
      {
        const std::size_t cell_index = it->attribute("cell_index").as_uint();
        const std::size_t local_entity
          = it->attribute("local_entity").as_uint();
        const std::string str_value = it->attribute("value").as_string();
        const std::size_t pos = str_value.find(',');
        dolfin_assert(pos != std::string::npos);
        const double re_value = std::stod(str_value.substr(0, pos));
        const double im_value
          = std::stod(str_value.substr(pos + 1, std::string::npos));
        const std::complex<double> value(re_value, im_value);
        mesh_value_collection.set_value(cell_index, local_entity, value);
      }
    }
    else
    {
      dolfin_error("XMLValueCollection.h",
                   "read mesh value collection from XML file",
                   "Unhandled value type \"%s\"", type.c_str());
    }
  }
}
