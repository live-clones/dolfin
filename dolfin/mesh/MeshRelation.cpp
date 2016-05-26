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
void MeshRelation::init(unsigned int dim, std::vector<std::size_t> map_back_to_src)
{
  dolfin_assert(src_mesh);
  dolfin_assert(dest_mesh);
  dolfin_assert(dim <= src_mesh->topology().dim());
  if (dest_to_src.size() < dim + 1)
  {
    dest_to_src.resize(dim + 1);
    dest_to_src_offset.resize(dim + 1);
  }

  auto& conn = dest_to_src[dim];

  if (!conn.empty())
  {
    dolfin_error("MeshRelation.cpp",
		 "initialise relation",
		 "Already initialised");
  }

  dolfin_assert(map_back_to_src.size() == dest_mesh->size(dim));
  conn.assign(map_back_to_src.begin(), map_back_to_src.end());
}
//-----------------------------------------------------------------------------
void MeshRelation::init(unsigned int dim, std::vector<std::vector<std::size_t>> map_back_to_src)
{
  dolfin_assert(src_mesh);
  dolfin_assert(dest_mesh);
  dolfin_assert(dim <= src_mesh->topology().dim());
  if (dest_to_src.size() < dim + 1)
  {
    dest_to_src.resize(dim + 1);
    dest_to_src_offset.resize(dim + 1);
  }

  auto& conn = dest_to_src[dim];
  auto& offset = dest_to_src_offset[dim];

  if (!conn.empty())
  {
    dolfin_error("MeshRelation.cpp",
		 "initialise relation",
		 "Already initialised");
  }

  dolfin_assert(map_back_to_src.size() == dest_mesh->size(dim));

  offset.push_back(0);
  for (unsigned int i = 0; i != map_back_to_src.size(); ++i)
  {
    conn.insert(conn.end(),
		map_back_to_src[i].begin(),
		map_back_to_src[i].end());
    offset.push_back(conn.size());
  }
}
//-----------------------------------------------------------------------------
unsigned int MeshRelation::num_entities(unsigned int dim, unsigned int index)
{
  dolfin_assert(dest_to_src_offset.size() > dim);
  auto& offset = dest_to_src_offset[dim];

  if (offset.empty())
    return 1;

  dolfin_assert(index < offset.size());
  return (offset[index + 1] - offset[index]);
}
//-----------------------------------------------------------------------------
unsigned int* MeshRelation::entities(unsigned int dim, unsigned int index)
{
  dolfin_assert(dest_to_src_offset.size() > dim);
  auto& offset = dest_to_src_offset[dim];
  auto& entity_indices = dest_to_src[dim];

  if (offset.empty())
  {
    dolfin_assert(index < dest_to_src.size());
    return &entity_indices[index];
  }

  dolfin_assert(index < offset.size());
  return &entity_indices[offset[index]];
}


