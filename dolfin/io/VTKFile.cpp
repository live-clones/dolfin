// Copyright (C) 2005-2009 Garth N. Wells
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
// Modified by Anders Logg 2005-2011
// Modified by Kristian Oelgaard 2006
// Modified by Martin Alnes 2008
// Modified by Niclas Jansson 2009
// Modified by Johannes Ring 2012
// Modified by Cian Wilson 2013
//
// First added:  2005-07-05
// Last changed: 2013-06-09

#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <boost/cstdint.hpp>
#include <boost/detail/endian.hpp>

#include "pugixml.hpp"

#include <dolfin/common/Timer.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include "Encoder.h"
#include "VTKWriter.h"
#include "VTKFile.h"

using namespace dolfin;

//----------------------------------------------------------------------------
VTKFile::VTKFile(const std::string filename, std::string encoding)
  : GenericFile(filename, "VTK"),
    _encoding(encoding), binary(false), compress(false)
{
  if (encoding != "ascii" && encoding != "base64" && encoding != "compressed")
  {
    dolfin_error("VTKFile.cpp",
                 "create VTK file",
                 "Unknown encoding (\"%s\"). "
                 "Known encodings are \"ascii\", \"base64\" and \"compressed\"",
                 encoding.c_str());
  }

  if (encoding == "ascii")
  {
    encode_string = "ascii";
    binary = false;
  }
  else if (encoding == "base64" || encoding == "compressed")
  {
    encode_string = "binary";
    binary = true;
    if (encoding == "compressed")
      compress = true;
  }
  else
  {
    dolfin_error("VTKFile.cpp",
                 "create VTK file",
                 "Unknown encoding (\"%s\"). "
                 "Known encodings are \"ascii\", \"base64\" and \"compressed\"",
                 encoding.c_str());
  }
}
//----------------------------------------------------------------------------
VTKFile::~VTKFile()
{
  // Do nothing
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const Mesh& mesh)
{
  write_mesh(mesh, counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const FunctionSpace& functionspace)
{
  write_functionspace(functionspace, counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const MeshFunction<bool>& meshfunction)
{
  mesh_function_write(meshfunction, counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const MeshFunction<std::size_t>& meshfunction)
{
  mesh_function_write(meshfunction, counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const MeshFunction<int>& meshfunction)
{
  mesh_function_write(meshfunction, counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const MeshFunction<double>& meshfunction)
{
  mesh_function_write(meshfunction, counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const Function& u)
{
  dolfin_assert(u.function_space()->mesh());
  const Mesh& mesh = *u.function_space()->mesh();

  std::vector<const GenericFunction*> us;
  us.push_back(&u);

  write(us, mesh, (double) counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::vector<const Function*>& us)
{
  std::vector<const Function*>::const_iterator u;
  std::vector<const GenericFunction*> usg;
  for (u = us.begin(); u != us.end(); u++)
  {
    usg.push_back(*u);
  }

  dolfin_assert(us[0]->function_space()->mesh());
  const Mesh& mesh = *us[0]->function_space()->mesh();

  write(usg, mesh, (double) counter);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const Mesh*, double> mesh)
{
  dolfin_assert(mesh.first);
  write_mesh(*(mesh.first), mesh.second);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const FunctionSpace*, double> functionspace)
{
  dolfin_assert(functionspace.first);
  write_functionspace(*(functionspace.first), functionspace.second);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const MeshFunction<int>*, double> f)
{
  dolfin_assert(f.first);
  mesh_function_write(*(f.first), f.second);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const MeshFunction<std::size_t>*,
                                         double> f)
{
  dolfin_assert(f.first);
  mesh_function_write(*(f.first), f.second);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const MeshFunction<double>*, double> f)
{
  dolfin_assert(f.first);
  mesh_function_write(*(f.first), f.second);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const MeshFunction<bool>*, double> f)
{
  dolfin_assert(f.first);
  mesh_function_write(*(f.first), f.second);
}
//----------------------------------------------------------------------------
void VTKFile::operator<<(const std::pair<const Function*, double> u)
{
  dolfin_assert(u.first);
  dolfin_assert(u.first->function_space()->mesh());
  const Mesh& mesh = *u.first->function_space()->mesh();

  std::vector<const GenericFunction*> us;
  us.push_back(u.first);

  write(us, mesh, u.second);
}
////----------------------------------------------------------------------------
//void VTKFile::operator<<(const std::pair<const std::vector<const Function*>, double> us)
//{
//  std::vector<const Function*>::const_iterator u;
//  std::vector<const GenericFunction*> usg;
//  for (u = us.first.begin(); u != us.first.end(); u++)
//  {
//    usg.push_back(*u);
//  }
//
//  dolfin_assert(us.first[0]->function_space()->mesh());
//  const Mesh& mesh = *us.first[0]->function_space()->mesh();
//
//  write(usg, mesh, us.second);
//}
//----------------------------------------------------------------------------
void VTKFile::write(const std::vector< std::shared_ptr<GenericFunction> >& us, const Mesh& mesh, double time)
{
  std::vector<const GenericFunction*> usp;
  std::vector<std::shared_ptr<GenericFunction> >::const_iterator u;
  for (u = us.begin(); u != us.end(); u++)
  {
    usp.push_back(&(**u));
  }
  write(usp, mesh, time);
}
//----------------------------------------------------------------------------
void VTKFile::write(const std::vector<const GenericFunction*>& us, const Mesh& mesh, double time)
{
  pugi::xml_document xml_doc;

  // Get MPI communicator
  const MPI_Comm mpi_comm = mesh.mpi_comm();

  // Get vtu file name and initialise
  std::string vtu_filename = init(xml_doc, mesh, mesh.topology().dim());

  // Write mesh
  VTKWriter::write_mesh(mesh, mesh.topology().dim(), xml_doc, binary,
                        compress);

  // Write results
  results_write(us, mesh, xml_doc);

  // Parallel-specific files
  const std::size_t num_processes = MPI::size(mpi_comm);
  if (num_processes > 1 && MPI::rank(mpi_comm) == 0)
  {
    std::string pvtu_filename = vtu_name(0, 0, counter, ".pvtu");
    pvtu_write(us, mesh, pvtu_filename);
    pvd_file_write(counter, time, pvtu_filename);
  }
  else if (num_processes == 1)
    pvd_file_write(counter, time, vtu_filename);

  // Finalise and write file
  finalize(xml_doc, vtu_filename);

  log(TRACE, "Saved functions to file %s in VTK format.", _filename.c_str());
}
//----------------------------------------------------------------------------
void VTKFile::write(const std::vector< std::shared_ptr<GenericFunction> >& us, const FunctionSpace& functionspace, double time)
{
  std::vector<const GenericFunction*> usp;
  std::vector<std::shared_ptr<GenericFunction> >::const_iterator u;
  for (u = us.begin(); u != us.end(); u++)
  {
    usp.push_back(&(**u));
  }
  write(usp, functionspace, time);
}
//----------------------------------------------------------------------------
void VTKFile::write(const std::vector<const GenericFunction*>& us, const FunctionSpace& functionspace, double time)
{
  dolfin_assert(functionspace.mesh());
  const Mesh& mesh = *functionspace.mesh();

  pugi::xml_document xml_doc;

  // Get MPI communicator
  const MPI_Comm mpi_comm = mesh.mpi_comm();

  // Get vtu file name and intialise
  std::string vtu_filename = init(xml_doc, functionspace, mesh.topology().dim());

  // Write mesh
  VTKWriter::write_mesh(functionspace, mesh.topology().dim(), xml_doc, binary,
                        compress);

  // Write results
  results_write(us, functionspace, xml_doc);

  // Parallel-specfic files
  const std::size_t num_processes = MPI::size(mpi_comm);
  if (num_processes > 1 && MPI::rank(mpi_comm) == 0)
  {
    std::string pvtu_filename = vtu_name(0, 0, counter, ".pvtu");
    pvtu_write(us, mesh, pvtu_filename);
    pvd_file_write(counter, time, pvtu_filename);
  }
  else if (num_processes == 1)
    pvd_file_write(counter, time, vtu_filename);

  // Finalise and write file
  finalize(xml_doc, vtu_filename);

  log(TRACE, "Saved functions to file %s in VTK format.", _filename.c_str());
}
//----------------------------------------------------------------------------
void VTKFile::write_mesh(const Mesh& mesh, double time)
{
  Timer t("Write mesh to PVD/VTK file");

  pugi::xml_document xml_doc;

  // Get MPI communicator
  const MPI_Comm mpi_comm = mesh.mpi_comm();

  // Get vtu file name and initialise out files
  std::string vtu_filename = init(xml_doc, mesh, mesh.topology().dim());

  // Write local mesh to vtu file
  VTKWriter::write_mesh(mesh, mesh.topology().dim(), xml_doc, binary, 
                        compress);

  // Parallel-specific files
  const std::size_t num_processes = MPI::size(mpi_comm);
  if (num_processes > 1 && MPI::rank(mpi_comm) == 0)
  {
    std::string pvtu_filename = vtu_name(0, 0, counter, ".pvtu");
    pvtu_write_mesh(pvtu_filename, num_processes);
    pvd_file_write(counter, time, pvtu_filename);
  }
  else if (num_processes == 1)
    pvd_file_write(counter, time, vtu_filename);

  // Finalise
  finalize(xml_doc, vtu_filename);

  log(TRACE, "Saved mesh %s (%s) to file %s in VTK format.",
      mesh.name().c_str(), mesh.label().c_str(), _filename.c_str());
}
//----------------------------------------------------------------------------
void VTKFile::write_functionspace(const FunctionSpace& functionspace, double time)
{
  dolfin_assert(functionspace.mesh());
  const Mesh& mesh = *functionspace.mesh();

  pugi::xml_document xml_doc;

  // Get MPI communicator
  const MPI_Comm mpi_comm = mesh.mpi_comm();

  // Get vtu file name and intialise out files
  std::string vtu_filename = init(xml_doc, functionspace, mesh.topology().dim());

  // Write local mesh to vtu file
  VTKWriter::write_mesh(functionspace, mesh.topology().dim(), xml_doc, binary, compress);

  // Parallel-specfic files
  const std::size_t num_processes = MPI::size(mpi_comm);
  if (num_processes > 1 && MPI::rank(mpi_comm) == 0)
  {
    std::string pvtu_filename = vtu_name(0, 0, counter, ".pvtu");
    pvtu_write_mesh(pvtu_filename, num_processes);
    pvd_file_write(counter, time, pvtu_filename);
  }
  else if (num_processes == 1)
    pvd_file_write(counter, time, vtu_filename);

  // Finalise
  finalize(xml_doc, vtu_filename);

  log(TRACE, "Saved functionspace on mesh %s (%s) to file %s in VTK format.",
      mesh.name().c_str(), mesh.label().c_str(), _filename.c_str());
}
//----------------------------------------------------------------------------
std::string VTKFile::init(pugi::xml_document& xml_doc, const Mesh& mesh, std::size_t cell_dim) const
{
  // Get MPI communicators
  const MPI_Comm mpi_comm = mesh.mpi_comm();

  // Get vtu file name
  std::string vtu_filename = vtu_name(MPI::rank(mpi_comm),
                                      MPI::size(mpi_comm),
                                      counter, ".vtu");

  // Number of cells and vertices
  const std::size_t num_cells = mesh.topology().ghost_offset(cell_dim);
  const std::size_t num_vertices = mesh.topology().ghost_offset(0);

  // Write headers
  vtk_header_open(num_vertices, num_cells, xml_doc);

  return vtu_filename;
}
//----------------------------------------------------------------------------
std::string VTKFile::init(pugi::xml_document& xml_doc, const FunctionSpace& functionspace, std::size_t cell_dim) const
{
  dolfin_assert(functionspace.element()->num_sub_elements() == 0);
  dolfin_assert(functionspace.element()->value_rank() == 0);

  std::shared_ptr<const GenericDofMap> dofmap = functionspace.dofmap();
  dolfin_assert(dofmap);

  const Mesh& mesh = *functionspace.mesh();

  // Get MPI communicators
  const MPI_Comm mpi_comm = mesh.mpi_comm();

  // Get vtu file name
  std::string vtu_filename = vtu_name(MPI::rank(mpi_comm),
                                      MPI::size(mpi_comm),
                                      counter, ".vtu");

  // Number of cells and vertices
  const std::size_t num_cells = mesh.topology().ghost_offset(cell_dim);

  // Create vector to hold dofs
  std::vector<la_index> dofs;
  dofs.reserve(num_cells*dofmap->max_element_dofs());

  ArrayView<const dolfin::la_index> cell_dofs;
  for (dolfin::CellIterator cell(mesh); !cell.end(); ++cell)       // loop over the cells in the mesh
  {
    // Check that cell is not a ghost
    dolfin_assert(!cell->is_ghost());

    cell_dofs = dofmap->cell_dofs((*cell).index());

    dofs.insert(dofs.end(), cell_dofs.begin(), cell_dofs.end());
  }

  // Sort dofs (required to later remove duplicates)
  std::sort(dofs.begin(), dofs.end());

  // Remove duplicates
  dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());

  // Write headers
  vtk_header_open(dofs.size(), num_cells, xml_doc);

  return vtu_filename;
}
//----------------------------------------------------------------------------
void VTKFile::finalize(pugi::xml_document& xml_doc, std::string vtu_filename)
{
  // Save file
  xml_doc.save_file(vtu_filename.c_str(), "  ");

  // Increase the number of times we have saved the object
  counter++;
}
//----------------------------------------------------------------------------
void VTKFile::results_write(const std::vector<const GenericFunction*>& us, const Mesh& mesh, pugi::xml_document& xml_doc) const
{

  std::vector<std::size_t> cell_counter(3,0);
  std::vector<std::size_t> point_counter(3,0);

  std::vector<const GenericFunction*>::const_iterator u;
  for (u = us.begin(); u != us.end(); u++)
  {
    // Get rank of Function
    const std::size_t rank = (*u)->value_rank();
    if (rank > 2)
    {
      dolfin_error("VTKFile.cpp",
                   "write data to VTK file",
                   "Only scalar, vector and tensor functions can be saved in VTK format");
    }

    // Get number of components
    const std::size_t dim = (*u)->value_size();

    // Check that function type can be handled
    if (rank == 1)
    {
      if (!(dim == 1 || dim == 2 || dim == 3))
      {
        dolfin_error("VTKFile.cpp",
                     "write data to VTK file",
                     "Don't know how to handle vector function with dimension other than 1, 2 or 3");
      }
    }
    else if (rank == 2)
    {
      if (!(dim == 4 || dim == 9))
      {
        dolfin_error("VTKFile.cpp",
                     "write data to VTK file",
                     "Don't know how to handle tensor function with dimension other than 4 or 9");
      }
    }


    const Function *uf = dynamic_cast<const Function*>(*u);
    if (uf)
    {
      std::size_t cell_based_dim = 1;
      dolfin_assert(uf->function_space()->mesh());
      dolfin_assert(uf->function_space()->dofmap());
      for (std::size_t i = 0; i < rank; i++)
        cell_based_dim *= uf->function_space()->mesh()->topology().dim();
      if ((uf->function_space()->dofmap()->max_element_dofs() == cell_based_dim) &&  // check if data is cell based
          (uf->function_space()->mesh()->num_cells()==mesh.num_cells()) && // also check it has the right number of cells
          (uf->function_space()->mesh()->topology().dim()==mesh.topology().dim())) // and the topo dim is the same
        VTKWriter::write_cell_data(*uf, xml_doc, binary, compress, cell_counter);
      else
        write_point_data(*uf, mesh, xml_doc, point_counter);
    }
    else
      write_point_data(**u, mesh, xml_doc, point_counter);
  }
}
//----------------------------------------------------------------------------
void VTKFile::results_write(const std::vector<const GenericFunction*>& us, const FunctionSpace& functionspace, pugi::xml_document& xml_doc) const
{

  dolfin_assert(functionspace.mesh());
  const Mesh& mesh = *functionspace.mesh();

  std::vector<std::size_t> cell_counter(3,0);
  std::vector<std::size_t> point_counter(3,0);

  std::vector<const GenericFunction*>::const_iterator u;
  for (u = us.begin(); u != us.end(); u++)
  {
    // Get rank of Function
    const std::size_t rank = (*u)->value_rank();
    if (rank > 2)
    {
      dolfin_error("VTKFile.cpp",
                   "write data to VTK file",
                   "Only scalar, vector and tensor functions can be saved in VTK format");
    }

    // Get number of components
    const std::size_t dim = (*u)->value_size();

    // Check that function type can be handled
    if (rank == 1)
    {
      if (!(dim == 1 || dim == 2 || dim == 3))
      {
        dolfin_error("VTKFile.cpp",
                     "write data to VTK file",
                     "Don't know how to handle vector function with dimension other than 1, 2 or 3");
      }
    }
    else if (rank == 2)
    {
      if (!(dim == 4 || dim == 9))
      {
        dolfin_error("VTKFile.cpp",
                     "write data to VTK file",
                     "Don't know how to handle tensor function with dimension other than 4 or 9");
      }
    }

    const Function *uf = dynamic_cast<const Function*>(*u);
    if (uf)
    {
      std::size_t cell_based_dim = 1;
      dolfin_assert(uf->function_space()->mesh());
      dolfin_assert(uf->function_space()->dofmap());
      for (std::size_t i = 0; i < rank; i++)
        cell_based_dim *= uf->function_space()->mesh()->topology().dim();
      if ((uf->function_space()->dofmap()->max_element_dofs() == cell_based_dim) &&  // check if data is cell based
          (uf->function_space()->mesh()->num_cells()==functionspace.mesh()->num_cells()) && // also check it has the right number of cells
          (uf->function_space()->mesh()->topology().dim()==functionspace.mesh()->topology().dim())) // and the topo dim is the same
        VTKWriter::write_cell_data(*uf, xml_doc, binary, compress, cell_counter);
      else
        write_point_data(*uf, functionspace, xml_doc, point_counter);
    }
    else
      write_point_data(**u, functionspace, xml_doc, point_counter);
  }
}
//----------------------------------------------------------------------------
void VTKFile::write_point_data(const GenericFunction& u, const Mesh& mesh,
                               pugi::xml_document& xml_doc, std::vector<std::size_t>& counter) const
{
  const std::size_t rank = u.value_rank();
  const std::size_t num_vertices = mesh.num_vertices();

  // Get number of components
  const std::size_t dim = u.value_size();

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  pugi::xml_node point_node, data_node, value_node;
  std::stringstream value;

  if ((counter[0]==0) && (counter[1]==0) && (counter[2]==0))
    point_node = piece_node.append_child("PointData");
  else
    point_node = piece_node.child("PointData");

  // Allocate memory for function values at vertices
  const std::size_t size = num_vertices*dim;
  std::vector<double> values(size);

  // Get function values at vertices
  u.compute_vertex_values(values, mesh);
  dolfin_assert(values.size() == size);

  if (rank == 0)
  {
    if (counter[rank]==0)
      point_node.append_attribute("Scalars") = u.name().c_str();

    data_node = point_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("format") = encode_string.c_str();

    counter[rank]++;
  }
  else if (rank == 1)
  {
    if (counter[rank]==0)
      point_node.append_attribute("Vectors") = u.name().c_str();

    data_node = point_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("NumberOfComponents") = "3";
    data_node.append_attribute("format") = encode_string.c_str();

    counter[rank]++;
  }
  else if (rank == 2)
  {
    if (counter[rank]==0)
      point_node.append_attribute("Tensors") = u.name().c_str();

    data_node = point_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("NumberOfComponents") = "9";
    data_node.append_attribute("format") = encode_string.c_str();

    counter[rank]++;
  }

  value.str("");
  if (_encoding == "ascii")
  {
    value << std::scientific;
    value << std::setprecision(16);
    for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    {
      if (rank == 1 && dim < 3)
      {
        // Append 0.0 to 2D vectors to make them 3D
        for(std::size_t i = 0; i < dim; i++)
          value << values[vertex->index() + i*num_vertices] << " ";
        for(std::size_t i = dim; i < 3; i++)
          value << 0.0 << "  ";
      }
      else if (rank == 2 && dim == 4)
      {
        // Pad 2D tensors with 0.0 to make them 3D
        for(std::size_t i = 0; i < 2; i++)
        {
          value << values[vertex->index() + (2*i + 0)*num_vertices] << " ";
          value << values[vertex->index() + (2*i + 1)*num_vertices] << " ";
          value << 0.0 << " ";
        }
        value << 0.0 << " ";
        value << 0.0 << " ";
        value << 0.0 << "  ";
      }
      else
      {
        // Write all components
        for(std::size_t i = 0; i < dim; i++)
          value << values[vertex->index() + i*num_vertices] << " ";
        value << " ";
      }
    }

    // Send to file
    value_node = data_node.append_child(pugi::node_pcdata);
    value_node.set_value(value.str().c_str());
  }
  else if (_encoding == "base64" || _encoding == "compressed")
  {
    // Number of zero paddings per point
    std::size_t padding_per_point = 0;
    std::vector<std::size_t> indicies(dim, 0);
    std::iota(indicies.begin(), indicies.end(), 0);
    if (rank == 1 && dim < 3)
      padding_per_point = 3-dim;
    else if (rank == 2 && dim == 4)
    {
      padding_per_point = 5;
      indicies[2] = 3;
      indicies[3] = 4;
    }

    // Number of data entries per point and total number
    const std::size_t num_data_per_point = dim + padding_per_point;
    const std::size_t num_total_data_points = num_vertices*num_data_per_point;

    std::vector<double> data(num_total_data_points, 0);
    for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    {
      const std::size_t index = vertex->index();
      for(std::size_t i = 0; i < dim; i++)
        data[index*num_data_per_point + indicies[i]] = values[index + i*num_vertices];
    }

    // Create encoded stream
    value << VTKWriter::encode_stream(data, compress) << std::endl;
    value_node = data_node.append_child(pugi::node_pcdata);
    value_node.set_value(value.str().c_str());
  }

}
//----------------------------------------------------------------------------
void VTKFile::write_point_data(const GenericFunction& u, const FunctionSpace& functionspace,
                               pugi::xml_document& xml_doc, std::vector<std::size_t>& counter) const
{
  const std::size_t rank = u.value_rank();

  // Get number of components
  const std::size_t dim = u.value_size();

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  pugi::xml_node point_node, data_node, value_node;
  std::stringstream value;

  if ((counter[0]==0) && (counter[1]==0) && (counter[2]==0))
    point_node = piece_node.append_child("PointData");
  else
    point_node = piece_node.child("PointData");

  std::shared_ptr<const GenericDofMap> dofmap = functionspace.dofmap();
  dolfin_assert(dofmap);

  const Mesh& mesh = *functionspace.mesh();

  // Create vector to hold dofs
  std::vector<la_index> dofs;
  dofs.reserve(mesh.num_cells()*dofmap->max_element_dofs());

  ArrayView<const dolfin::la_index> cell_dofs;
  for (dolfin::CellIterator cell(mesh); !cell.end(); ++cell)       // loop over the cells in the mesh
  {
    cell_dofs = dofmap->cell_dofs((*cell).index());

    dofs.insert(dofs.end(), cell_dofs.begin(), cell_dofs.end());
  }

  // Sort dofs (required to later remove duplicates)
  std::sort(dofs.begin(), dofs.end());

  // Remove duplicates
  dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
  
  // Allocate memory for function values at functionspace dofs
  Function uf(std::make_shared<FunctionSpace>(functionspace));
  std::vector< std::vector<double> > values(dim);
  for (std::vector< std::vector<double> >::iterator v = values.begin();
                                                    v != values.end();
                                                    v++)
  {
    v->resize(dofs.size());
  }

  if (rank == 0)
  {
    if (counter[rank]==0)
      point_node.append_attribute("Scalars") = u.name().c_str();

    data_node = point_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("format") = encode_string.c_str();

    uf.interpolate(u);
    uf.vector()->get_local(values[0].data(), dofs.size(), dofs.data());
    std::vector<double>::iterator it;
    for (it = values[0].begin(); it != values[0].end(); ++it)
    {
      if (std::abs(*it) < DOLFIN_EPS)
        *it = 0.0;
    }

    counter[rank]++;
  }
  else if (rank == 1)
  {
    if (counter[rank]==0)
      point_node.append_attribute("Vectors") = u.name().c_str();

    data_node = point_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("NumberOfComponents") = "3";
    data_node.append_attribute("format") = encode_string.c_str();

    const Function* ul = dynamic_cast<const Function*>(&u);
    if(ul && (ul->function_space()->element()->num_sub_elements()==dim))
    {
      for (std::size_t i = 0; i < dim; i++)
      {
        uf.interpolate((*ul)[i]);
        uf.vector()->get_local(values[i].data(), dofs.size(), dofs.data());
        std::vector<double>::iterator it;
        for (it = values[i].begin(); it != values[i].end(); ++it)
        {
          if (std::abs(*it) < DOLFIN_EPS)
            *it = 0.0;
        }
      }
    }
    else // FIXME: Functions with no sub elements are very costly to evaluate here
    {
      for (std::size_t i = 0; i < dim; i++)
      {
        VTKExpression ue(i, &u);
        uf.interpolate(ue);
        uf.vector()->get_local(values[i].data(), dofs.size(), dofs.data());
        std::vector<double>::iterator it;
        for (it = values[i].begin(); it != values[i].end(); ++it)
        {
          if (std::abs(*it) < DOLFIN_EPS)
            *it = 0.0;
        }
      }
    }

    counter[rank]++;
  }
  else if (rank == 2)
  {
    if (counter[rank]==0)
      point_node.append_attribute("Tensors") = u.name().c_str();

    data_node = point_node.append_child("DataArray");
    data_node.append_attribute("type") = "Float64";
    data_node.append_attribute("Name") = u.name().c_str();
    data_node.append_attribute("NumberOfComponents") = "9";
    data_node.append_attribute("format") = encode_string.c_str();

    const Function* ul = dynamic_cast<const Function*>(&u);
    if(ul && (ul->function_space()->element()->num_sub_elements() > 0))
    {
      const bool symmetric = (ul->function_space()->element()->num_sub_elements() != ul->value_size());
      const std::size_t dim0 = ul->value_dimension(0);
      const std::size_t dim1 = ul->value_dimension(1);
      for (std::size_t i = 0; i < dim0; i++)
      {
        for (std::size_t j = 0; j < dim1; j++)
        {
          std::size_t k = i*dim1 + j;
          std::size_t kf = k;
          if (symmetric)
            if (j >= i)
              kf = i*dim1 + j - (i*(i+1))/2;
            else
              kf = j*dim1 + i - (j*(j+1))/2;
          uf.interpolate((*ul)[kf]);
          uf.vector()->get_local(values[k].data(), dofs.size(), dofs.data());
          std::vector<double>::iterator it;
          for (it = values[k].begin(); it != values[k].end(); ++it)
          {
            if (std::abs(*it) < DOLFIN_EPS)
              *it = 0.0;
          }
        }
      }
    }
    else // FIXME: Functions with no sub elements are very costly to evaluate here
    {
      const std::size_t dim0 = u.value_dimension(0);
      const std::size_t dim1 = u.value_dimension(1);
      for (std::size_t i = 0; i < dim0; i++)
      {
        for (std::size_t j = 0; j < dim1; j++)
        {
          std::size_t k = i*dim1 + j;
          VTKExpression ue(k, &u);
          uf.interpolate(ue);
          uf.vector()->get_local(values[k].data(), dofs.size(), dofs.data());
          std::vector<double>::iterator it;
          for (it = values[k].begin(); it != values[k].end(); ++it)
          {
            if (std::abs(*it) < DOLFIN_EPS)
              *it = 0.0;
          }
        }
      }
    }

    counter[rank]++;
  }

  std::size_t lsize = values[0].size();

  value.str("");
  if (_encoding == "ascii")
  {
    value << std::scientific;
    value << std::setprecision(16);
    for (std::size_t index = 0; index < lsize; index++)
    {
      if (rank == 1 && dim < 3)
      {
        // Append 0.0 to 2D vectors to make them 3D
        for(std::size_t i = 0; i < dim; i++)
          value << values[i][index] << " ";
        for(std::size_t i = dim; i < 3; i++)
          value << 0.0 << "  ";
      }
      else if (rank == 2 && dim == 4)
      {
        // Pad 2D tensors with 0.0 to make them 3D
        for(std::size_t i = 0; i < 2; i++)
        {
          value << values[2*i][index] << " ";
          value << values[2*i+1][index] << " ";
          value << 0.0 << " ";
        }
        value << 0.0 << " ";
        value << 0.0 << " ";
        value << 0.0 << "  ";
      }
      else
      {
        // Write all components
        for(std::size_t i = 0; i < dim; i++)
          value << values[i][index] << " ";
        value << " ";
      }
    }

    // Send to file
    value_node = data_node.append_child(pugi::node_pcdata);
    value_node.set_value(value.str().c_str());
  }
  else if (_encoding == "base64" || _encoding == "compressed")
  {
    // Number of zero paddings per point
    std::size_t padding_per_point = 0;
    std::vector<std::size_t> indicies(dim, 0);
    std::iota(indicies.begin(), indicies.end(), 0);
    if (rank == 1 && dim < 3)
      padding_per_point = 3-dim;
    else if (rank == 2 && dim == 4)
    {
      padding_per_point = 5;
      indicies[2] = 3;
      indicies[3] = 4;
    }

    // Number of data entries per point and total number
    const std::size_t num_data_per_point = dim + padding_per_point;
    const std::size_t num_total_data_points = lsize*num_data_per_point;

    std::vector<double> data(num_total_data_points, 0);
    for (std::size_t index = 0; index < lsize; index++)
    {
      for(std::size_t i = 0; i < dim; i++)
        data[index*num_data_per_point + indicies[i]] = values[i][index];
    }

    // Create encoded stream
    value << VTKWriter::encode_stream(data, compress) << std::endl;
    value_node = data_node.append_child(pugi::node_pcdata);
    value_node.set_value(value.str().c_str());
  }

}
//----------------------------------------------------------------------------
void VTKFile::pvd_file_write(std::size_t step, double time,
                             std::string fname)
{
  pugi::xml_document xml_doc;
  if (step == 0)
  {
    pugi::xml_node vtk_node = xml_doc.append_child("VTKFile");
    vtk_node.append_attribute("type") = "Collection";
    vtk_node.append_attribute("version") = "0.1";
    vtk_node.append_child("Collection");
  }
  else
  {
    pugi::xml_parse_result result = xml_doc.load_file(_filename.c_str());
    if (!result)
    {
      dolfin_error("VTKFile.cpp",
                   "write data to VTK file",
                   "XML parsing error when reading from existing file");
    }
  }

  // Remove directory path from name for pvd file
  const std::string fname_strip = strip_path(fname);

  // Get Collection node
  pugi::xml_node xml_collections = xml_doc.child("VTKFile").child("Collection");
  dolfin_assert(xml_collections);

  // Append data set
  pugi::xml_node dataset_node = xml_collections.append_child("DataSet");
  dataset_node.append_attribute("timestep") = time;
  dataset_node.append_attribute("part") = "0";
  dataset_node.append_attribute("file") = fname_strip.c_str();

  // Save file
  xml_doc.save_file(_filename.c_str(), "  ");
}
//----------------------------------------------------------------------------
void VTKFile::pvtu_write_mesh(pugi::xml_node xml_node) const
{
  // Vertex data
  pugi::xml_node vertex_data_node = xml_node.append_child("PPoints");
  pugi::xml_node data_node = vertex_data_node.append_child("PDataArray");
  data_node.append_attribute("type") = "Float64";
  data_node.append_attribute("NumberOfComponents") = "3";

  // Cell data
  pugi::xml_node cell_data_node = xml_node.append_child("PCellData");

  data_node = cell_data_node.append_child("PDataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "connectivity";

  data_node = cell_data_node.append_child("PDataArray");
  data_node.append_attribute("type") = "UInt32";
  data_node.append_attribute("Name") = "offsets";

  data_node = cell_data_node.append_child("PDataArray");
  data_node.append_attribute("type") = "UInt8";
  data_node.append_attribute("Name") = "types";
}
//----------------------------------------------------------------------------
void VTKFile::pvtu_write_function(std::size_t dim, std::size_t rank,
                                  const std::string data_location,
                                  const std::string name,
                                  const std::string fname,
                                  std::size_t num_processes,
                                  std::vector<std::size_t>& point_counter,
                                  std::vector<std::size_t>& cell_counter) const
{
  // Create xml doc
  pugi::xml_document xml_doc;
  pugi::xml_node vtk_node, grid_node, data_node;

  if ((cell_counter[0]==0) && (cell_counter[1]==0) && (cell_counter[2]==0) &&
      (point_counter[0]==0) && (point_counter[1]==0) && (point_counter[2]==0))
  {
    vtk_node = xml_doc.append_child("VTKFile");
    vtk_node.append_attribute("type") = "PUnstructuredGrid";
    vtk_node.append_attribute("version") = "0.1";
    grid_node = vtk_node.append_child("PUnstructuredGrid");
    grid_node.append_attribute("GhostLevel") = 0;

    // Mesh
    pvtu_write_mesh(grid_node);

    // Write vtu file list
    for(std::size_t i = 0; i < num_processes; i++)
    {
      const std::string tmp_string = strip_path(vtu_name(i, num_processes, counter, ".vtu"));
      pugi::xml_node piece_node = grid_node.append_child("Piece");
      piece_node.append_attribute("Source") = tmp_string.c_str();
    }

  }
  else
  {
    pugi::xml_parse_result result = xml_doc.load_file(fname.c_str());
    if (!result)
    {
      dolfin_error("VTKFile.cpp",
                   "write data to pvtu file",
                   "XML parsing error when reading from existing file");
    }
    grid_node = xml_doc.child("VTKFile").child("PUnstructuredGrid");
    dolfin_assert(grid_node);
  }

  // Add function data
  if (data_location == "point")
  {
    if ((point_counter[0]==0) && (point_counter[1]==0) && (point_counter[2]==0))
      data_node = grid_node.append_child("PPointData");
    else
      data_node = grid_node.child("PPointData");
  }
  else if (data_location == "cell")
  {
    if ((cell_counter[0]==0) && (cell_counter[1]==0) && (cell_counter[2]==0))
      data_node = grid_node.append_child("PCellData");
    else
      data_node = grid_node.child("PCellData");
  }

  // Get type based on rank
  std::size_t num_components = 0;
  if (rank == 0)
  {
    if (data_location == "point")
    {
      if (point_counter[rank] == 0)
        data_node.append_attribute("Scalars") = name.c_str();
      point_counter[rank]++;
    }
    else if (data_location == "cell") 
    {
      if (cell_counter[rank] == 0)
        data_node.append_attribute("Scalars") = name.c_str();
      cell_counter[rank]++;
    }

    num_components = 0;
  }
  else if (rank == 1)
  {
    if (!(dim == 2 || dim == 3))
    {
      dolfin_error("VTKFile.cpp",
                   "write data to VTK file",
                   "Don't know how to handle vector function with dimension other than 2 or 3");
    }

    if (data_location == "point")
    {
      if (point_counter[rank] == 0)
        data_node.append_attribute("Vectors") = name.c_str();
      point_counter[rank]++;
    }
    else if (data_location == "cell") 
    {
      if (cell_counter[rank] == 0)
        data_node.append_attribute("Vectors") = name.c_str();
      cell_counter[rank]++;
    }

    num_components = 3;
  }
  else if (rank == 2)
  {
    if (!(dim == 4 || dim == 9))
    {
      dolfin_error("VTKFile.cpp",
                   "write data to VTK file",
                   "Don't know how to handle tensor function with dimension other than 4 or 9");
    }

    if (data_location == "point")
    {
      if (point_counter[rank] == 0)
        data_node.append_attribute("Tensors") = name.c_str();
      point_counter[rank]++;
    }
    else if (data_location == "cell") 
    {
      if (cell_counter[rank] == 0)
        data_node.append_attribute("Tensors") = name.c_str();
      cell_counter[rank]++;
    }

    num_components = 9;
  }
  else
  {
    dolfin_error("VTKFile.cpp",
                 "write data to VTK file",
                 "Cannot handle XML output of rank %d", rank);
  }

  pugi::xml_node data_array_node = data_node.append_child("PDataArray");
  data_array_node.append_attribute("type") = "Float64";
  data_array_node.append_attribute("Name") = name.c_str();
  if (num_components>0)
    data_array_node.append_attribute("NumberOfComponents")
      = (unsigned int) num_components;

  xml_doc.save_file(fname.c_str(), "  ");
}
//----------------------------------------------------------------------------
void VTKFile::pvtu_write_mesh(const std::string fname,
                              const std::size_t num_processes) const
{
  // Create xml doc
  pugi::xml_document xml_doc;
  pugi::xml_node vtk_node = xml_doc.append_child("VTKFile");
  vtk_node.append_attribute("type") = "PUnstructuredGrid";
  vtk_node.append_attribute("version") = "0.1";
  pugi::xml_node grid_node = vtk_node.append_child("PUnstructuredGrid");
  grid_node.append_attribute("GhostLevel") = 0;

  // Mesh
  pvtu_write_mesh(grid_node);

  // Write vtu file list
  for (std::size_t i = 0; i < num_processes; i++)
  {
    const std::string tmp_string = strip_path(vtu_name(i, num_processes,
                                                       counter, ".vtu"));
    pugi::xml_node piece_node = grid_node.append_child("Piece");
    piece_node.append_attribute("Source") = tmp_string.c_str();
  }

  xml_doc.save_file(fname.c_str(), "  ");
}
//----------------------------------------------------------------------------
void VTKFile::pvtu_write(const std::vector<const GenericFunction*>& us, const Mesh& mesh, const std::string fname) const
{
  std::vector<std::size_t> cell_counter(3,0);
  std::vector<std::size_t> point_counter(3,0);

  const std::size_t num_processes = MPI::size(mesh.mpi_comm());

  std::vector<const GenericFunction*>::const_iterator u;
  for (u = us.begin(); u != us.end(); u++)
  {
    const std::size_t rank = (*u)->value_rank();
    if (rank > 2)
    {
      dolfin_error("VTKFile.cpp",
                   "write data to VTK file",
                   "Only scalar, vector and tensor functions can be saved in VTK format");
    }

    // Get number of components
    const std::size_t dim = (*u)->value_size();

    // Test for cell-based element type
    std::string data_type = "point";

    const Function* uf = dynamic_cast<const Function*>(*u);
    if (uf)
    {
      std::size_t cell_based_dim = 1;
      dolfin_assert(uf->function_space()->mesh());
      dolfin_assert(uf->function_space()->dofmap());
      for (std::size_t i = 0; i < rank; i++)
        cell_based_dim *= uf->function_space()->mesh()->topology().dim();
      if ((uf->function_space()->dofmap()->max_element_dofs() == cell_based_dim) &&  // check if data is cell based
          (uf->function_space()->mesh()->num_cells()==mesh.num_cells()) && // also check it has the right number of cells
          (uf->function_space()->mesh()->topology().dim()==mesh.topology().dim())) // and the topo dim is the same
        data_type = "cell";
    }

    pvtu_write_function(dim, rank, data_type, (*u)->name(), fname, num_processes,
                        point_counter, cell_counter);
  }

}
//----------------------------------------------------------------------------
void VTKFile::vtk_header_open(std::size_t num_points, std::size_t num_cells,
                              pugi::xml_document& xml_doc) const
{
  // Write headers
  pugi::xml_node vtk_node = xml_doc.append_child("VTKFile");
  vtk_node.append_attribute("type") = "UnstructuredGrid";
  vtk_node.append_attribute("version") = "0.1";
  if (encode_string == "binary")
  {
    #if defined BOOST_LITTLE_ENDIAN
    vtk_node.append_attribute("byte_order") = "LittleEndian";
    #elif defined BOOST_BIG_ENDIAN
    vtk_node.append_attribute("byte_order") = "BigEndian";
    #else
    dolfin_error("VTKFile.cpp",
                 "write data to VTK file",
                 "Unable to determine the endianness of the machine for VTK binary output");
    #endif
  }

  // Compression string
  if (_encoding == "compressed")
    vtk_node.append_attribute("compressor") = "vtkZLibDataCompressor";
  
  pugi::xml_node grid_node = vtk_node.append_child("UnstructuredGrid");

  pugi::xml_node piece_node = grid_node.append_child("Piece");
  piece_node.append_attribute("NumberOfPoints") = (unsigned int) num_points;
  piece_node.append_attribute("NumberOfCells") = (unsigned int) num_cells;

}
//----------------------------------------------------------------------------
std::string VTKFile::vtu_name(const int process, const int num_processes,
                              const int counter, std::string ext) const
{
  std::string filestart, extension;
  std::ostringstream fileid, newfilename;

  fileid.fill('0');
  fileid.width(6);

  filestart.assign(_filename, 0, _filename.find_last_of("."));
  extension.assign(_filename, _filename.find_last_of("."), _filename.size());

  fileid << counter;

  // Add process number to .vtu file name
  std::string proc = "";
  if (num_processes > 1)
  {
    std::ostringstream _p;
    _p << "_p" << process << "_";
    proc = _p.str();
  }
  newfilename << filestart << proc << fileid.str() << ext;

  return newfilename.str();
}
//----------------------------------------------------------------------------
template<typename T>
void VTKFile::mesh_function_write(T& meshfunction, double time)
{
  const Mesh& mesh = *meshfunction.mesh();
  const std::size_t cell_dim = meshfunction.dim();

  pugi::xml_document xml_doc;

  // Update vtu file name and clear file
  std::string vtu_filename = init(xml_doc, mesh, cell_dim);

  // Write mesh
  VTKWriter::write_mesh(mesh, cell_dim, xml_doc, binary, compress);

  pugi::xml_node piece_node = xml_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

  pugi::xml_node cell_node, data_node, value_node;
  std::stringstream value;

  cell_node = piece_node.append_child("CellData");
  cell_node.append_attribute("Scalars") = meshfunction.name().c_str();

  data_node = cell_node.append_child("DataArray");
  data_node.append_attribute("type") = "Float64";
  data_node.append_attribute("Name") = meshfunction.name().c_str();
  data_node.append_attribute("format") = "ascii";
  // Open file
  value.str("");
  // Write data
  for (MeshEntityIterator cell(mesh, cell_dim); !cell.end(); ++cell)
    value << meshfunction[cell->index()] << " ";

  value_node = data_node.append_child(pugi::node_pcdata);
  value_node.set_value(value.str().c_str());

  // Parallel-specific files
  const std::size_t num_processes = MPI::size(mesh.mpi_comm());
  const std::size_t process_number = MPI::rank(mesh.mpi_comm());
  if (num_processes > 1 && process_number == 0)
  {
    std::string pvtu_filename = vtu_name(0, 0, counter, ".pvtu");
    std::vector<std::size_t> cell_counter(3,0), point_counter(3,0);
    pvtu_write_function(1, 0, "cell", meshfunction.name(), pvtu_filename,
                        num_processes, point_counter, cell_counter);
    pvd_file_write(counter, time, pvtu_filename);
  }
  else if (num_processes == 1)
    pvd_file_write(counter, time, vtu_filename);

  // Finalise
  finalize(xml_doc, vtu_filename);

  log(TRACE, "Saved mesh function %s (%s) to file %s in VTK format.",
      mesh.name().c_str(), mesh.label().c_str(), _filename.c_str());
}
//----------------------------------------------------------------------------
std::string VTKFile::strip_path(std::string file) const
{
  std::string fname;
  fname.assign(file, _filename.find_last_of("/") + 1, file.size());
  return fname;
}
//----------------------------------------------------------------------------
