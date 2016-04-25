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
// Modified by Anders Logg 2006.
// Modified by Niclas Jansson 2009.
//
// First added:  2005-07-05
// Last changed: 2013-03-11

#ifndef __VTK_FILE_H
#define __VTK_FILE_H

#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include "GenericFile.h"
#include "pugixml.hpp"
#include <boost/shared_ptr.hpp>
#include <dolfin/function/Expression.h>

namespace pugi
{
  class xml_node;
}

namespace dolfin
{

  /// This class supports the output of meshes and functions in VTK
  /// XML format for visualisation purposes. It is not suitable to
  /// checkpointing as it may decimate some data.

  class VTKFile : public GenericFile
  {
  public:

    VTKFile(const std::string filename, std::string encoding);
    ~VTKFile();

    void operator<< (const Mesh& mesh);
    void operator<< (const FunctionSpace& functionspace);
    void operator<< (const MeshFunction<bool>& meshfunction);
    void operator<< (const MeshFunction<std::size_t>& meshfunction);
    void operator<< (const MeshFunction<int>& meshfunction);
    void operator<< (const MeshFunction<double>& meshfunction);
    void operator<< (const Function& u);
    void operator<< (const std::vector<const Function*>& us);
    void operator<< (const std::pair<const Mesh*, double> mesh);
    void operator<< (const std::pair<const FunctionSpace*, double> functionspace);
    void operator<< (const std::pair<const MeshFunction<int>*, double> f);
    void
      operator<< (const std::pair<const MeshFunction<std::size_t>*, double> f);
    void operator<< (const std::pair<const MeshFunction<double>*, double> f);
    void operator<< (const std::pair<const MeshFunction<bool>*, double> f);
    void operator<< (const std::pair<const Function*, double> u);
    //void operator<< (const std::pair<const std::vector<const Function*>, double> us);

    void write(const std::vector<std::shared_ptr<GenericFunction> >& us, const Mesh& mesh, double time);
    void write(const std::vector<const GenericFunction*>& us, const Mesh& mesh, double time);
    void write(const std::vector<std::shared_ptr<GenericFunction> >& us, const FunctionSpace& functionspace, double time);
    void write(const std::vector<const GenericFunction*>& us, const FunctionSpace& functionspace, double time);

  protected:

    void write_mesh(const Mesh& mesh, double time);

    void write_functionspace(const FunctionSpace& functionspace, double time);

    std::string init(pugi::xml_document& xml_doc, const Mesh& mesh, std::size_t dim) const;

    std::string init(pugi::xml_document& xml_doc, const FunctionSpace& functionspace, std::size_t dim) const;

    void finalize(pugi::xml_document& xml_doc, std::string vtu_filename);

    void results_write(const std::vector<const GenericFunction*>& us, const Mesh& mesh, pugi::xml_document& xml_doc) const;

    void results_write(const std::vector<const GenericFunction*>& us, const FunctionSpace& functionspace, pugi::xml_document& xml_doc) const;

    void write_point_data(const GenericFunction& u, const Mesh& mesh,
                          pugi::xml_document& xml_doc, std::vector<std::size_t>& counter) const;

    void write_point_data(const GenericFunction& u, const FunctionSpace& functionspace,
                          pugi::xml_document& xml_doc, std::vector<std::size_t>& counter) const;

    void pvd_file_write(std::size_t step, double time, std::string fname);


    void pvtu_write_function(std::size_t dim, std::size_t rank,
                             const std::string data_location,
                             const std::string name,
                             const std::string fname,
                             std::size_t num_processes,
                             std::vector<std::size_t>& point_counter,
                             std::vector<std::size_t>& cell_counter) const;

    void pvtu_write_mesh(const std::string pvtu_filename,
                         const std::size_t num_processes) const;

    void pvtu_write(const std::vector<const GenericFunction*>& us, const Mesh& mesh, const std::string pvtu_filename) const;

    void vtk_header_open(std::size_t num_points, std::size_t num_cells, 
                         pugi::xml_document& xml_doc) const;

    void vtk_header_close(std::string file) const;

    std::string vtu_name(const int process, const int num_processes,
                         const int counter, std::string ext) const;

    template<typename T>
    void mesh_function_write(T& meshfunction, double time);

    // Strip path from file
    std::string strip_path(std::string file) const;

  private:

    void pvtu_write_mesh(pugi::xml_node xml_node) const;

    // File encoding
    const std::string _encoding;
    std::string encode_string;

    bool binary;
    bool compress;

  };

  class VTKExpression : public Expression
  {

  public:

    VTKExpression(const std::size_t& dim, const GenericFunction* u) : _dim(dim), _u(u)
    {
      dolfin_assert(_dim < _u->value_size());
    }

    void eval(Array<double>& values, const Array<double>& x, const ufc::cell &cell) const
    {
      Array<double> _u_values(_u->value_size());
      _u->eval(_u_values, x, cell);
      values[0] = _u_values[_dim];
    }

  private:

    const std::size_t _dim;

    const GenericFunction* _u;

  };

}

#endif
