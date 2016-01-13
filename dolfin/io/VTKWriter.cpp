// Copyright (C) 2010 Garth N. Wells
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
// Modified by Anders Logg 2011
// Modified by Johannes Ring 2012
//
// First added:  2010-07-19
// Last changed: 2012-09-14

#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <boost/detail/endian.hpp>
#include "pugixml.hpp"

#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include "Encoder.h"
#include "VTKWriter.h"

using namespace dolfin;

//----------------------------------------------------------------------------
void VTKWriter::write_mesh(const Mesh& mesh, std::size_t cell_dim,
                           pugi::xml_document& xml_doc, bool binary, bool compress)
{
  if (binary)
    write_base64_mesh(mesh, cell_dim, xml_doc, compress);
  else
    write_ascii_mesh(mesh, cell_dim, xml_doc);
}
//----------------------------------------------------------------------------
void VTKWriter::write_mesh(const FunctionSpace& functionspace, std::size_t cell_dim,
                           pugi::xml_document& xml_doc, bool binary, bool compress)
{
  if (binary)
    write_base64_mesh(functionspace, cell_dim, xml_doc, compress);
  else
    write_ascii_mesh(functionspace, cell_dim, xml_doc);
}
//----------------------------------------------------------------------------
void VTKWriter::write_cell_data(const Function& u, pugi::xml_document& xml_doc,
                                bool binary, bool compress, std::vector<std::size_t>& counter)
{
  // For brevity
  dolfin_assert(u.function_space()->mesh());
  dolfin_assert(u.function_space()->dofmap());
  const Mesh& mesh = *u.function_space()->mesh();
  const GenericDofMap& dofmap = *u.function_space()->dofmap();
  const std::size_t tdim = mesh.topology().dim();
  const std::size_t num_cells = mesh.topology().ghost_offset(tdim);

  std::string encode_string;
  if (!binary)
    encode_string = "ascii";
  else
    encode_string = "binary";

  // Get rank of Function
  const std::size_t rank = u.value_rank();
  if(rank > 2)
  {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Don't know how to handle functions with rank greater than 2");
  }

  // Get number of components
  const std::size_t data_dim = u.value_size();

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  pugi::xml_node cell_node, data_node, value_node;
  std::stringstream value;

  if ((counter[0]==0) && (counter[1]==0) && (counter[2]==0))
    cell_node = piece_node.append_child("CellData");
  else
    cell_node = piece_node.child("CellData");

  // Write headers
  if (rank == 0)
  {
    if (counter[rank]==0)
      cell_node.append_attribute("Scalars") = u.name().c_str();

    data_node = cell_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("format") = encode_string.c_str();

    counter[rank]++;
  }
  else if (rank == 1)
  {
    if(!(data_dim == 2 || data_dim == 3))
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Don't know how to handle vector function with dimension other than 2 or 3");
    }

    if (counter[rank]==0)
      cell_node.append_attribute("Vectors") = u.name().c_str();

    data_node = cell_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("NumberOfComponents") = "3";
    data_node.append_attribute("format") = encode_string.c_str();

    counter[rank]++;
  }
  else if (rank == 2)
  {
    if(!(data_dim == 4 || data_dim == 9))
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Don't know how to handle tensor function with dimension other than 4 or 9");
    }

    if (counter[rank]==0)
      cell_node.append_attribute("Tensors") = u.name().c_str();

    data_node = cell_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("NumberOfComponents") = "9";
    data_node.append_attribute("format") = encode_string.c_str();

    counter[rank]++;
  }

  // Allocate memory for function values at cell centres
  const std::size_t size = num_cells*data_dim;

  // Build lists of dofs and create map
  std::vector<dolfin::la_index> dof_set;
  std::vector<std::size_t> offset(size + 1);
  std::vector<std::size_t>::iterator cell_offset = offset.begin();
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Tabulate dofs
    const ArrayView<const dolfin::la_index>
      dofs = dofmap.cell_dofs(cell->index());
    for(std::size_t i = 0; i < dofmap.num_element_dofs(cell->index()); ++i)
      dof_set.push_back(dofs[i]);

    // Add local dimension to cell offset and increment
    *(cell_offset + 1)
      = *(cell_offset) + dofmap.num_element_dofs(cell->index());
    ++cell_offset;
  }

  // Get  values
  std::vector<double> values(dof_set.size());
  dolfin_assert(u.vector());
  u.vector()->get_local(values.data(), dof_set.size(), dof_set.data());

  // Get cell data
  value.str("");
  if (!binary)
    value << ascii_cell_data(mesh, offset, values, data_dim, rank);
  else
  {
    value << base64_cell_data(mesh, offset, values, data_dim, rank, compress) 
          << std::endl;
  }
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());
  
}
//----------------------------------------------------------------------------
std::string VTKWriter::ascii_cell_data(const Mesh& mesh,
                                       const std::vector<std::size_t>& offset,
                                       const std::vector<double>& values,
                                       std::size_t data_dim, std::size_t rank)
{
  std::ostringstream ss;
  ss << std::scientific;
  ss << std::setprecision(16);
  std::vector<std::size_t>::const_iterator cell_offset = offset.begin();
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    if (rank == 1 && data_dim == 2)
    {
      // Append 0.0 to 2D vectors to make them 3D
      ss << values[*cell_offset] << "  " << values[*cell_offset + 1] << " "
         << 0.0;
    }
    else if (rank == 2 && data_dim == 4)
    {
      // Pad with 0.0 to 2D tensors to make them 3D
      for(std::size_t i = 0; i < 2; i++)
      {
        ss << values[*cell_offset + 2*i] << " ";
        ss << values[*cell_offset + 2*i + 1] << " ";
        ss << 0.0 << " ";
      }
      ss << 0.0 << " ";
      ss << 0.0 << " ";
      ss << 0.0;
    }
    else
    {
      // Write all components
      for (std::size_t i = 0; i < data_dim; i++)
        ss << values[*cell_offset + i] << " ";
    }
    ss << "  ";
    ++cell_offset;
  }

  return ss.str();
}
//----------------------------------------------------------------------------
std::string VTKWriter::base64_cell_data(const Mesh& mesh,
                                        const std::vector<std::size_t>& offset,
                                        const std::vector<double>& values,
                                        std::size_t data_dim, std::size_t rank,
                                        bool compress)
{
  const std::size_t num_cells = mesh.num_cells();

  // Number of zero paddings per point
  std::size_t padding_per_point = 0;
  if (rank == 1 && data_dim == 2)
    padding_per_point = 1;
  else if (rank == 2 && data_dim == 4)
    padding_per_point = 5;

  // Number of data entries per point and total number
  const std::size_t num_data_per_point = data_dim + padding_per_point;
  const std::size_t num_total_data_points = num_cells*num_data_per_point;

  std::vector<std::size_t>::const_iterator cell_offset = offset.begin();
  std::vector<double> data(num_total_data_points, 0);
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    const std::size_t index = cell->index();
    for(std::size_t i = 0; i < data_dim; i++)
      data[index*num_data_per_point + i] = values[*cell_offset + i];
    ++cell_offset;
  }

  return encode_stream(data, compress);
}
//----------------------------------------------------------------------------
void VTKWriter::write_ascii_mesh(const Mesh& mesh, std::size_t cell_dim,
                                 pugi::xml_document& xml_doc)
{
  const std::size_t num_cells = mesh.topology().ghost_offset(cell_dim);
  const std::size_t num_cell_vertices = mesh.type().num_vertices(cell_dim);

  // Get VTK cell type
  const std::size_t _vtk_cell_type = vtk_cell_type(mesh, cell_dim);

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  pugi::xml_node data_node, value_node;
  std::stringstream value;

  // Write vertex positions
  pugi::xml_node points_node = piece_node.append_child("Points");

  data_node = points_node.append_child("DataArray");
  data_node.append_attribute("type") = "Float64";
  data_node.append_attribute("NumberOfComponents") = "3";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Point p = v->point();
    value << p.x() << " " << p.y() << " " <<  p.z() << "  ";
  }
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write cell connectivity
  pugi::xml_node cells_node = piece_node.append_child("Cells");
  
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "connectivity";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  std::unique_ptr<CellType>
    celltype(CellType::create(mesh.type().entity_type(cell_dim)));
  const std::vector<unsigned int> perm = celltype->vtk_mapping();
  for (MeshEntityIterator c(mesh, cell_dim); !c.end(); ++c)
  {
    for (unsigned int i = 0; i != c->num_entities(0); ++i)
      value << c->entities(0)[perm[i]] << " ";
    value << " ";
  }
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write offset into connectivity array for the end of each cell
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "offsets";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (std::size_t offsets = 1; offsets <= num_cells; offsets++)
    value << offsets*num_cell_vertices << " ";
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write cell type
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt8";
  data_node.append_attribute("Name") = "types";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (std::size_t types = 0; types < num_cells; types++)
    value << _vtk_cell_type << " ";
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

}
//----------------------------------------------------------------------------
void VTKWriter::write_ascii_mesh(const FunctionSpace& functionspace, std::size_t cell_dim,
                                 pugi::xml_document& xml_doc)
{
  // Assuming we have a non-mixed element
  dolfin_assert(functionspace.element()->num_sub_elements() == 0);
  dolfin_assert(functionspace.element()->value_rank() == 0);

  const Mesh& mesh = *functionspace.mesh();
  dolfin_assert(functionspace.dofmap());
  const GenericDofMap& dofmap = *functionspace.dofmap();
  dolfin_assert(functionspace.element());
  std::shared_ptr<const dolfin::FiniteElement> element = functionspace.element();

  const std::size_t num_cells = mesh.topology().size(cell_dim);
  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t cdim = dofmap.max_element_dofs();

  // Get VTK cell type
  const std::size_t _vtk_cell_type = vtk_cell_type(functionspace, cell_dim);

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  // Coordinates for dofs
  boost::multi_array<double, 2> cellcoords(boost::extents[cdim][gdim]);
  std::vector<double> vertexcoords;
  std::vector<double> dofcoords(3, 0.0);
  std::unordered_map<dolfin::la_index, std::vector<double> > coordinatemap;
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    c->get_coordinate_dofs(vertexcoords);
    element->tabulate_dof_coordinates(cellcoords, vertexcoords, *c);
    const ArrayView<const dolfin::la_index> dofs = dofmap.cell_dofs(c->index());
    // Loop over all dofs on cell
    for (std::size_t i = 0; i < dofmap.num_element_dofs(c->index()); ++i)
    {
      for (std::size_t j = 0; j < gdim; ++j)
        dofcoords[j] = cellcoords[i][j];
      coordinatemap.insert( std::unordered_map<dolfin::la_index, std::vector<double> >::value_type(dofs[i], dofcoords) );
    }
  }

  std::map<dolfin::la_index, std::vector<double> > ordered_coordinatemap(coordinatemap.begin(), coordinatemap.end());

  std::map<dolfin::la_index, dolfin::la_index> local_dofmap;
  dolfin::la_index index = 0;
  for (std::map<dolfin::la_index, std::vector<double> >::const_iterator c = ordered_coordinatemap.begin(); 
         c != ordered_coordinatemap.end(); c++)
  {
    local_dofmap.insert( std::map<dolfin::la_index, dolfin::la_index>::value_type(c->first, index++) );
  }

  pugi::xml_node data_node, value_node;
  std::stringstream value;

  // Write dof positions
  pugi::xml_node points_node = piece_node.append_child("Points");

  data_node = points_node.append_child("DataArray");
  data_node.append_attribute("type") = "Float64";
  data_node.append_attribute("NumberOfComponents") = "3";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (std::map<dolfin::la_index, std::vector<double> >::const_iterator c = ordered_coordinatemap.begin(); 
         c != ordered_coordinatemap.end(); c++)
  {
    value << c->second[0] << " " << c->second[1] << " " << c->second[2] << "  ";
  }
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  const std::vector<std::size_t> _vtk_cell_order = vtk_cell_order(functionspace, cell_dim);
  dolfin_assert(_vtk_cell_order.size()>=cdim);

  // Write cell connectivity
  pugi::xml_node cells_node = piece_node.append_child("Cells");
  
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "connectivity";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    const ArrayView<const dolfin::la_index> dofs = dofmap.cell_dofs(c->index());
    for (std::size_t i = 0; i < dofmap.num_element_dofs(c->index()); ++i)
      value << local_dofmap[dofs[_vtk_cell_order[i]]] << " ";
    value << " ";
  }
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write offset into connectivity array for the end of each cell
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "offsets";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (std::size_t offsets = 1; offsets <= num_cells; offsets++)
    value << offsets*(dofmap.num_element_dofs(offsets-1)) << " ";
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write cell type
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt8";
  data_node.append_attribute("Name") = "types";
  data_node.append_attribute("format") = "ascii";

  value.str("");
  for (std::size_t types = 0; types < num_cells; types++)
    value << _vtk_cell_type << " ";
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

}
//-----------------------------------------------------------------------------
void VTKWriter::write_base64_mesh(const Mesh& mesh, std::size_t cell_dim,
                                  pugi::xml_document& xml_doc, bool compress)
{
  const std::size_t num_cells = mesh.topology().size(cell_dim);
  const std::size_t num_cell_vertices = mesh.type().num_vertices(cell_dim);

  // Get VTK cell type
  const boost::uint8_t _vtk_cell_type = vtk_cell_type(mesh, cell_dim);

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  pugi::xml_node data_node, value_node;
  std::stringstream value;

  // Write vertex positions
  pugi::xml_node points_node = piece_node.append_child("Points");

  data_node = points_node.append_child("DataArray");
  data_node.append_attribute("type") = "Float64";
  data_node.append_attribute("NumberOfComponents") = "3";
  data_node.append_attribute("format") = "binary";

  std::vector<double> vertex_data(3*mesh.num_vertices());
  std::vector<double>::iterator vertex_entry = vertex_data.begin();
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    const Point p = v->point();
    *vertex_entry++ = p.x();
    *vertex_entry++ = p.y();
    *vertex_entry++ = p.z();
  }
  // Create encoded stream
  value.str("");
  value <<  encode_stream(vertex_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write cell connectivity
  pugi::xml_node cells_node = piece_node.append_child("Cells");
  
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "connectivity";
  data_node.append_attribute("format") = "binary";

  const int size = num_cells*num_cell_vertices;
  std::vector<boost::uint32_t> cell_data(size);
  std::vector<boost::uint32_t>::iterator cell_entry = cell_data.begin();

  std::unique_ptr<CellType>
    celltype(CellType::create(mesh.type().entity_type(cell_dim)));
  const std::vector<unsigned int> perm = celltype->vtk_mapping();
  for (MeshEntityIterator c(mesh, cell_dim); !c.end(); ++c)
  {
    for (unsigned int i = 0; i != c->num_entities(0); ++i)
      *cell_entry++ = c->entities(0)[perm[i]];
  }

  // Create encoded stream
  value.str("");
  value << encode_stream(cell_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write offset into connectivity array for the end of each cell
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "offsets";
  data_node.append_attribute("format") = "binary";

  std::vector<boost::uint32_t> offset_data(num_cells);
  std::vector<boost::uint32_t>::iterator offset_entry = offset_data.begin();
  for (std::size_t offsets = 1; offsets <= num_cells; offsets++)
    *offset_entry++ = offsets*num_cell_vertices;

  // Create encoded stream
  value.str("");
  value << encode_stream(offset_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write cell type
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt8";
  data_node.append_attribute("Name") = "types";
  data_node.append_attribute("format") = "binary";

  std::vector<boost::uint8_t> type_data(num_cells);
  std::vector<boost::uint8_t>::iterator type_entry = type_data.begin();
  for (std::size_t types = 0; types < num_cells; types++)
    *type_entry++ = _vtk_cell_type;

  // Create encoded stream
  value.str("");
  value << encode_stream(type_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

}
//-----------------------------------------------------------------------------
void VTKWriter::write_base64_mesh(const FunctionSpace& functionspace, std::size_t cell_dim,
                                  pugi::xml_document& xml_doc, bool compress)
{
  // Assuming we have a non-mixed element
  dolfin_assert(functionspace.element()->num_sub_elements() == 0);
  dolfin_assert(functionspace.element()->value_rank() == 0);

  const Mesh& mesh = *functionspace.mesh();
  dolfin_assert(functionspace.dofmap());
  const GenericDofMap& dofmap = *functionspace.dofmap();
  dolfin_assert(functionspace.element());
  std::shared_ptr<const dolfin::FiniteElement> element = functionspace.element();

  const std::size_t num_cells = mesh.topology().size(cell_dim);
  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t cdim = dofmap.max_element_dofs();

  // Get VTK cell type
  const boost::uint8_t _vtk_cell_type = vtk_cell_type(functionspace, cell_dim);

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  // Coordinates for dofs
  boost::multi_array<double, 2> cellcoords(boost::extents[cdim][gdim]);
  std::vector<double> vertexcoords;
  std::vector<double> dofcoords(3, 0.0);
  std::unordered_map<dolfin::la_index, std::vector<double> > coordinatemap;
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    c->get_coordinate_dofs(vertexcoords);
    element->tabulate_dof_coordinates(cellcoords, vertexcoords, *c);
    const ArrayView<const dolfin::la_index> dofs = dofmap.cell_dofs(c->index());
    // Loop over all dofs on cell
    for (std::size_t i = 0; i < dofmap.num_element_dofs(c->index()); ++i)
    {
      for (std::size_t j = 0; j < gdim; ++j)
        dofcoords[j] = cellcoords[i][j];
      coordinatemap.insert( std::unordered_map<dolfin::la_index, std::vector<double> >::value_type(dofs[i], dofcoords) );
    }
  }

  std::map<dolfin::la_index, std::vector<double> > ordered_coordinatemap(coordinatemap.begin(), coordinatemap.end());

  std::map<dolfin::la_index, dolfin::la_index> local_dofmap;
  dolfin::la_index index = 0;
  for (std::map<dolfin::la_index, std::vector<double> >::const_iterator c = ordered_coordinatemap.begin(); 
         c != ordered_coordinatemap.end(); c++)
  {
    local_dofmap.insert( std::map<dolfin::la_index, dolfin::la_index>::value_type(c->first, index++) );
  }

  pugi::xml_node data_node, value_node;
  std::stringstream value;

  // Write dof positions
  pugi::xml_node points_node = piece_node.append_child("Points");

  data_node = points_node.append_child("DataArray");
  data_node.append_attribute("type") = "Float64";
  data_node.append_attribute("NumberOfComponents") = "3";
  data_node.append_attribute("format") = "binary";

  std::vector<double> vertex_data(3*ordered_coordinatemap.size());
  std::vector<double>::iterator vertex_entry = vertex_data.begin();
  for (std::map<dolfin::la_index, std::vector<double> >::const_iterator c = ordered_coordinatemap.begin();
        c != ordered_coordinatemap.end(); c++)
  {
    *vertex_entry++ = c->second[0];
    *vertex_entry++ = c->second[1];
    *vertex_entry++ = c->second[2];
  }
  // Create encoded stream
  value.str("");
  value <<  encode_stream(vertex_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  const std::vector<std::size_t> _vtk_cell_order = vtk_cell_order(functionspace, cell_dim);
  dolfin_assert(_vtk_cell_order.size()>=cdim);

  // Write cell connectivity
  pugi::xml_node cells_node = piece_node.append_child("Cells");
  
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "connectivity";
  data_node.append_attribute("format") = "binary";

  std::vector<boost::uint32_t> cell_data(cdim*num_cells);
  std::vector<boost::uint32_t>::iterator cell_entry = cell_data.begin();
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    const ArrayView<const dolfin::la_index> dofs = dofmap.cell_dofs(c->index());
    for (std::size_t i = 0; i < dofmap.num_element_dofs(c->index()); ++i)
      *cell_entry++ = local_dofmap[dofs[_vtk_cell_order[i]]];
  }

  // Create encoded stream
  value.str("");
  value << encode_stream(cell_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write offset into connectivity array for the end of each cell
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "offsets";
  data_node.append_attribute("format") = "binary";

  std::vector<boost::uint32_t> offset_data(num_cells);
  std::vector<boost::uint32_t>::iterator offset_entry = offset_data.begin();
  for (std::size_t offsets = 1; offsets <= num_cells; offsets++)
    *offset_entry++ = offsets*(dofmap.num_element_dofs(offsets-1));

  // Create encoded stream
  value.str("");
  value << encode_stream(offset_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Write cell type
  data_node = cells_node.append_child("DataArray");
  data_node.append_attribute("type") = "UInt8";
  data_node.append_attribute("Name") = "types";
  data_node.append_attribute("format") = "binary";

  std::vector<boost::uint8_t> type_data(num_cells);
  std::vector<boost::uint8_t>::iterator type_entry = type_data.begin();
  for (std::size_t types = 0; types < num_cells; types++)
    *type_entry++ = _vtk_cell_type;

  // Create encoded stream
  value.str("");
  value << encode_stream(type_data, compress) << std::endl;
  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

}
//----------------------------------------------------------------------------
boost::uint8_t VTKWriter::vtk_cell_type(const Mesh& mesh,
                                        std::size_t cell_dim)
{
  // Get cell type
  CellType::Type cell_type = mesh.type().entity_type(cell_dim);

  // Determine VTK cell type
  boost::uint8_t vtk_cell_type = 0;
  if (cell_type == CellType::tetrahedron)
    vtk_cell_type = 10;
  else if (cell_type == CellType::hexahedron)
    vtk_cell_type = 12;
  else if (cell_type == CellType::quadrilateral)
    vtk_cell_type = 9;
  else if (cell_type == CellType::triangle)
    vtk_cell_type = 5;
  else if (cell_type == CellType::interval)
    vtk_cell_type = 3;
  else if (cell_type == CellType::point)
    vtk_cell_type = 1;
  else
  {
    dolfin_error("VTKWriter.cpp",
                 "write data to VTK file",
                 "Unknown cell type (%d)", cell_type);
  }

  return vtk_cell_type;
}
//----------------------------------------------------------------------------
boost::uint8_t VTKWriter::vtk_cell_type(const FunctionSpace& functionspace, std::size_t cell_dim)
{
  // Make sure this function space is based on a Lagrange element
  // FIXME: better way of checking this?
  dolfin_assert(functionspace.element()->signature().find("Lagrange")!=std::string::npos);

  const Mesh& mesh = *functionspace.mesh();
  // Get cell type
  CellType::Type cell_type = mesh.type().cell_type();
  if (mesh.topology().dim() == cell_dim)
    cell_type = mesh.type().cell_type();
  else if (mesh.topology().dim() - 1 == cell_dim)
    cell_type = mesh.type().facet_type();
  else if (cell_dim == 1)
    cell_type = CellType::interval;
  else if (cell_dim == 0)
    cell_type = CellType::point;
  else
  {
    dolfin_error("VTKWriter.cpp",
                 "write data to VTK file",
                 "Can only handle cells, cell facets or points with VTK output for now");
  }

  const std::size_t sdim = functionspace.element()->space_dimension();
  // Determine VTK cell type
  boost::uint8_t vtk_cell_type = 0;
  if (cell_type == CellType::tetrahedron)
  {
    if (sdim==10)
      vtk_cell_type = 24;
    else if (sdim==4)
      vtk_cell_type = 10;
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for tetradedron", sdim);
    }
  }
  else if (cell_type == CellType::triangle)
  {
    if (sdim==6)
      vtk_cell_type = 22;
    else if (sdim==3)
      vtk_cell_type = 5;
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for triangle", sdim);
    }
  }
  else if (cell_type == CellType::interval)
  {
    if (sdim==3)
      vtk_cell_type = 21;
    else if (sdim==2)
      vtk_cell_type = 3;
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for interval", sdim);
    }
  }
  else if (cell_type == CellType::point)
    if (sdim==1)
    {
      vtk_cell_type = 1;
    }
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for point", sdim);
    }
  else
  {
    dolfin_error("VTKWriter.cpp",
                 "write data to VTK file",
                 "Unknown cell type (%d)", cell_type);
  }

  return vtk_cell_type;
}
//----------------------------------------------------------------------------
std::vector<std::size_t> VTKWriter::vtk_cell_order(const FunctionSpace& functionspace, std::size_t cell_dim)
{
  // Make sure this function space is based on a Lagrange element
  // FIXME: better way of checking this?
  dolfin_assert(functionspace.element()->signature().find("Lagrange")!=std::string::npos);

  const Mesh& mesh = *functionspace.mesh();
  // Get cell type
  CellType::Type cell_type = mesh.type().cell_type();
  if (mesh.topology().dim() == cell_dim)
    cell_type = mesh.type().cell_type();
  else if (mesh.topology().dim() - 1 == cell_dim)
    cell_type = mesh.type().facet_type();
  else
  {
    dolfin_error("VTKWriter.cpp",
                 "write data to VTK file",
                 "Can only handle cells and cell facets with VTK output for now");
  }

  const std::size_t sdim = functionspace.element()->space_dimension();
  // Determine VTK cell type
  std::vector<std::size_t> vtk_cell_order;
  if (cell_type == CellType::tetrahedron)
  {
    if (sdim==10)
    {
      std::size_t cell_order[] = {0,1,2,3,9,6,8,7,5,4};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else if (sdim==4)
    {
      std::size_t cell_order[] = {0,1,2,3};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for tetradedron", sdim);
    }
  }
  else if (cell_type == CellType::triangle)
  {
    if (sdim==6)
    {
      std::size_t cell_order[] = {0,1,2,5,3,4};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else if (sdim==3)
    {
      std::size_t cell_order[] = {0,1,2};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for triangle", sdim);
    }
  }
  else if (cell_type == CellType::interval)
  {
    if (sdim==3)
    {
      std::size_t cell_order[] = {0,1,2};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else if (sdim==2)
    {
      std::size_t cell_order[] = {0,1};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for interval", sdim);
    }
  }
  else if (cell_type == CellType::point)
    if (sdim==1)
    {
      std::size_t cell_order[] = {0};
      vtk_cell_order.assign(cell_order, cell_order+sdim);
    }
    else
    {
      dolfin_error("VTKWriter.cpp",
                   "write data to VTK file",
                   "Unsupported space dimension (%d) for point", sdim);
    }
  else
  {
    dolfin_error("VTKWriter.cpp",
                 "write data to VTK file",
                 "Unknown cell type (%d)", cell_type);
  }

  return vtk_cell_order;
}
//----------------------------------------------------------------------------
