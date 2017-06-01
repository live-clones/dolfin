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
// First added:  2013-05-30
// Last changed: 2013-05-30

#include <dolfin/mesh/Mesh.h>
#include "MeshEntityIntersection.h"
#include "intersect.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::shared_ptr<const MeshEntityIntersection>
dolfin::intersect(const Mesh& mesh, const Point& point)
{
  return intersect(mesh, point, mesh.geometry().dim());
}
//-----------------------------------------------------------------------------
std::shared_ptr<const MeshEntityIntersection>
dolfin::intersect(const Mesh& mesh, const Point& point,
                  const std::size_t t_dim)
{
  return std::shared_ptr<const MeshEntityIntersection>
      (new MeshEntityIntersection(mesh, point, t_dim));
}
//-----------------------------------------------------------------------------
std::shared_ptr<const MeshEntityIntersection>
dolfin::intersect(const Mesh& mesh, const Point& x1, const Point& x2)
{
  return intersect(mesh, x1, x2, mesh.geometry().dim());
}
//-----------------------------------------------------------------------------
std::shared_ptr<const MeshEntityIntersection>
dolfin::intersect(const Mesh& mesh, const Point& x1, const Point& x2,
                  const std::size_t t_dim)
{
  return std::shared_ptr<const MeshEntityIntersection>
      (new MeshEntityIntersection(mesh, x1, x2, t_dim));
}
//-----------------------------------------------------------------------------
