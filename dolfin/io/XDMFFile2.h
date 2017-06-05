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

#ifndef __DOLFIN_XDMFFILE2_H
#define __DOLFIN_XDMFFILE2_H

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#ifdef HAS_HDF5
#include <hdf5.h>
#else
typedef int hid_t;
#endif

#include <dolfin/common/MPI.h>
#include <dolfin/common/Variable.h>
#include "PugiXDMFXMLDocument.h"

namespace boost
{
namespace filesystem
{
class path;
}
}

namespace pugi
{
class xml_node;
class xml_document;
}

namespace dolfin
{

// Forward declarations
class Function;
#ifdef HAS_HDF5
class HDF5File;
#endif
class LocalMeshData;
class Mesh;
template<typename T>
class MeshFunction;
template<typename T>
class MeshValueCollection;
class Point;
class XDMFxml;

/// Read and write Mesh, Function, MeshFunction and other objects in XDMF

/// This class supports the output of meshes and functions in XDMF
/// (http://www.xdmf.org) format. It creates an XML file that
/// describes the data and points to a HDF5 file that stores the
/// actual problem data. Output of data in parallel is supported.
///
/// XDMF is not suitable for checkpointing as it may decimate some
/// data.

class XDMFFile2 : public Variable
{
public:

  /// File encoding type
  enum class Encoding
  {
    HDF5, ASCII
  };

  /// Constructor
  XDMFFile2(const std::string filename)
      : XDMFFile2(MPI_COMM_WORLD, filename)
  {}

  /// Constructor
  XDMFFile2(MPI_Comm comm, const std::string filename);

  /// Destructor
  ~XDMFFile2();

  /// Close the file
  ///
  /// This closes any open HDF5 files. In ASCII mode the XML file is
  /// closed each time it is written to or read from, so close() has
  /// no effect.
  ///
  /// From Python you can also use XDMFFile as a context manager:
  ///
  ///     with XDMFFile(mpi_comm_world(), 'name.xdmf') as xdmf:
  ///         xdmf.write(mesh)
  ///
  /// The file is automatically closed at the end of the with block
  void close();

  /// Save a mesh to XDMF format, either using an associated HDF5
  /// file, or storing the data inline as XML Create function on
  /// given function space
  ///
  /// @param    mesh (_Mesh_)
  ///         A mesh to save.
  /// @param    encoding (_Encoding_)
  ///         Encoding to use: HDF5 or ASCII
  ///
  void write(const Mesh &mesh, Encoding encoding = XDMFFile2::_default_encoding);

private:

  // MPI communicator
  MPI_Comm _mpi_comm;

  // HDF5 dependent variables
  // Default encoding = HDF5 iff it is supported
#ifdef HAS_HDF5
  std::unique_ptr<HDF5File> _hdf5_file;
  static const Encoding _default_encoding = Encoding::HDF5;
#else
  const bool _hdf5_file = false;
  static const Encoding _default_encoding = Encoding::ASCII;
#endif

  // Cached filename
  const std::string _filename;

  // Counter for time series
  std::size_t _counter;

  // The XML document currently representing the XDMF
  // which needs to be kept open for time series etc.
  std::unique_ptr<PugiXDMFXMLDocument> _xdmf_xml_doc;

  // Register default XDMF parameters into global dolfin pool
  void add_dolfin_parameters();

  // Dolfin to VTK cell type conversion
  // FIXME: probably should be in CellType
  std::string vtk_cell_type_str(CellType::Type cell_type, const size_t degree);

  // Check whether the requested encoding is supported
  void check_encoding(Encoding encoding) const;

  // Checks the mesh
  void check_mesh(const Mesh &mesh) const;

  std::set<unsigned int> compute_nonlocal_entities(const Mesh &mesh, int cell_dim);

};

}

#endif
