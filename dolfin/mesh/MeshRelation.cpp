// Copyright (C) 2016 Chris Richardson
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

#include<map>
#include "Mesh.h"

#include "MeshRelation.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshRelation::init(unsigned int dim, std::vector<std::size_t> map)
{
  dolfin_assert(src_mesh);
  dolfin_assert(dest_mesh);
  dolfin_assert(dim <= src_mesh->topology().dim());
  if (src_to_dest.size() < dim + 1)
    src_to_dest.resize(dim + 1);

  auto& conn = src_to_dest[dim];

  if (!conn.empty())
  {
    dolfin_error("MeshRelation.cpp",
		 "initialise relation",
		 "Already initialised");
  }

  dolfin_assert(map.size() == src_mesh->size(dim));
  conn.assign(map.begin(), map.end());
}
//-----------------------------------------------------------------------------
void MeshRelation::init(unsigned int dim, std::vector<std::vector<std::size_t>> map)
{
  dolfin_assert(src_mesh);
  dolfin_assert(dest_mesh);
  dolfin_assert(dim <= src_mesh->topology().dim());
  if (src_to_dest.size() < dim + 1)
    src_to_dest.resize(dim + 1);

  auto& conn = src_to_dest[dim];

  if (!conn.empty())
  {
    dolfin_error("MeshRelation.cpp",
		 "initialise relation",
		 "Already initialised");
  }

  dolfin_assert(map.size() == src_mesh->size(dim));

  for (unsigned int i = 0; i != map.size(); ++i)
    conn.insert(conn.end(), map[i].begin(), map[i].end());
}
//-----------------------------------------------------------------------------
