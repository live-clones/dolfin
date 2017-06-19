// Copyright (C) 2017 Nate Sime and Chris Richardson
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

#include <dolfin/function/Function.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include "Point.h"
#include "CollisionDetection.h"

#include "GeometricContact.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::vector<Point> GeometricContact::create_deformed_segment_volume_3d(Mesh& mesh, std::size_t facet_index, const Function& u)
{
  Facet facet(mesh, facet_index);
  Point X1 = Vertex(mesh, facet.entities(0)[0]).point();
  Point X2 = Vertex(mesh, facet.entities(0)[1]).point();
  Point X3 = Vertex(mesh, facet.entities(0)[2]).point();
  Point x1 = X1 + Point(u(X1));
  Point x2 = X2 + Point(u(X2));
  Point x3 = X3 + Point(u(X3));


  // Cast out prism defined by three tetrahedra from the master facet to its displaced configuration
  return {X1, X2, X3, x1,
          X2, X3, x1, x2,
          X3, x1, x2, x3};
}
//-----------------------------------------------------------------------------
bool GeometricContact::check_tet_set_collision(const std::vector<Point>& tet_set1, const std::vector<Point>& tet_set2)
{
  dolfin_assert(tet_set1.size() == tet_set2.size());
  dolfin_assert(tet_set1.size() %4 == 0);
  const std::size_t ntets = tet_set1.size() / 4;

  for (unsigned int i = 0; i < ntets; ++i)
    for (unsigned int j = 0; i < ntets; ++i)
    {
      if (CollisionDetection::collides_tetrahedron_tetrahedron(tet_set1[i*4],
                                                               tet_set1[i*4+1],
                                                               tet_set1[i*4+2],
                                                               tet_set1[i*4+3],
                                                               tet_set2[j*4],
                                                               tet_set2[j*4+1],
                                                               tet_set2[j*4+2],
                                                               tet_set2[j*4+3]))
        return true;
    }

  return false;
}
//-----------------------------------------------------------------------------
std::map<std::size_t, std::vector<std::size_t>>
  GeometricContact::contact_surface_map_volume_sweep_3d(Mesh& mesh, Function& u,
  const std::vector<std::size_t>& master_facets, const std::vector<std::size_t>& slave_facets)
{
  // Construct a dictionary mapping master facets to their collided slave counterparts.

  // This algorithm find which slave facets are contained in the volume swept out by
  // the master surface and its displacement.

  // :param mesh: The mesh
  // :param u: The (3D) solution vector field
  // :param master_facets: list of master facet indices
  // :param slave_facets: list of slave facet indices
  // :return: mapping master facet index to slave facet indices

  std::map<std::size_t, std::vector<std::size_t>> master_to_slave;

  for (auto &mf : master_facets)
  {
    std::vector<Point> master_tet_set = create_deformed_segment_volume_3d(mesh, mf, u);
    for (auto &sf : slave_facets)
    {
      std::vector<Point> slave_tet_set = create_deformed_segment_volume_3d(mesh, sf, u);
      if (check_tet_set_collision(master_tet_set, slave_tet_set))
        master_to_slave[mf].push_back(sf);
    }
  }

  return master_to_slave;
}
