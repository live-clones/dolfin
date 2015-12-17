// Copyright (C) 2015 Ettie Unwin
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


#include <cmath>
#include <dolfin/common/constants.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/math/basic.h>

#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include "Plane.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
/// Checks to see if an edge intersects a plane and returns T/F and the
/// intersection point.  If no intersection returns (0, 0, 0).
std::pair<bool, Point> Plane::intersection( const Edge& e )
{
  std::pair<bool, Point>  result = {false, Point(0.0, 0.0, 0.0)};

  const Mesh& mesh = e.mesh();
  Point p0 = Vertex(mesh, e.entities(0)[0]).point();
  Point p1 = Vertex(mesh, e.entities(0)[1]).point();

  if (std::signbit(p0.dot(normal()) - _d) !=
      std::signbit(p1.dot(normal()) - _d) )
  {
    result.first = true;
    // FIXME: what if p0.n == p1.n ?
    double lambda = (_d - p1.dot(_n))/(p0.dot(_n) - p1.dot(_n));
    result.second = lambda*p0 + (1-lambda)*p1;
  }
  else
    result.first = false;

  return result;
}
//-----------------------------------------------------------------------------
