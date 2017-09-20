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

#ifndef __JSON_FILE_H
#define __JSON_FILE_H

#include "GenericFile.h"

namespace dolfin
{

  /// This class implements output of tables to json format(JSON)

  class JSONFile : public GenericFile
  {
  public:

    /// Constructor
    explicit JSONFile(const std::string filename);

    /// Destructor
    ~JSONFile();

    /// Write table to file in a json format
    void write(const Table& table);

  };

}

#endif
