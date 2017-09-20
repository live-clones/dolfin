// Copyright (C) 2017 Jan Hybs
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

#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <dolfin/common/MPI.h>
#include <dolfin/log/log.h>
#include <dolfin/parameter/GlobalParameters.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include "JSONFile.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
JSONFile::JSONFile(const std::string filename): GenericFile(filename, "JSON")
{
  // Do nothing
}
//-----------------------------------------------------------------------------
JSONFile::~JSONFile()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void JSONFile::write(const Table& table)
{

  auto it_double(table.dvalues.cbegin());
  auto it_string(table.values.cbegin());

  const auto row_count = table.rows.size();
  const auto col_count = table.cols.size();

  boost::property_tree::ptree root, rows;

  root.put<std::string>("name", table.name().c_str());
  for (std::size_t i = 0; i < row_count; ++i)
  {
    boost::property_tree::ptree row, cols;
    row.put<std::string>("name", table.rows[i].c_str());

    for (std::size_t j = 0; j < col_count; ++j)
    {
      boost::property_tree::ptree col;
      col.put<std::string>("key", table.cols[j].c_str());

      const auto key = std::make_pair(table.rows[i], table.cols[j]);
      it_double = table.dvalues.find(key);
      if (it_double != table.dvalues.end())
      {
        col.put<std::string>("type", "double");
        col.put("value", it_double->second);
      }
      else
      {
        it_string = table.values.find(key);
        col.put<std::string>("type", "double");
        col.put<std::string>("value", it_string->second.c_str());
      }
      cols.push_back( std::make_pair("", col) );
    }
    row.put_child("cols", cols);
    rows.push_back(std::make_pair("", row));
  }
  root.put_child("rows", rows);

  // Open file
  std::ofstream fp(_filename.c_str(), std::ios_base::out);
  boost::property_tree::write_json(fp, root, true);
}
//-----------------------------------------------------------------------------
