// Copyright (C) 2013 Anders Logg
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
// First added:  2013-04-18
// Last changed: 2013-08-28

#include <dolfin/log/LogStream.h>
#include <dolfin/mesh/Cell.h>
#include "BoundingBoxTree.h"
#include "MeshEntityIntersection.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshEntityIntersection::MeshEntityIntersection(const Mesh& mesh,
                                               const Point& point,
                                               const std::size_t t_dim)
{
  // Build bounding box tree
  BoundingBoxTree tree;
  tree.build(mesh, t_dim);

  // Compute intersection
  _intersected_entities = tree.compute_entity_collisions(point);
}
//-----------------------------------------------------------------------------
MeshEntityIntersection::MeshEntityIntersection(const Mesh& mesh,
                                               const Point& x1,
                                               const Point& x2,
                                               const std::size_t t_dim)
{
  // Build bounding box tree
  BoundingBoxTree tree;
  tree.build(mesh, t_dim);

  // Compute intersection
  _intersected_entities = tree.compute_entity_collisions(x1, x2);
}
//-----------------------------------------------------------------------------
