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

#ifndef __GEOMETRIC_CONTACT_H
#define __GEOMETRIC_CONTACT_H

#include <map>
#include <vector>

namespace dolfin
{

  // Forward declarations
  class Mesh;
  class Point;
  class Facet;
  class Function;

  /// This class implements ...

  class GeometricContact
  {
  public:

    /// Return map from master facets to possible colliding slave facets
    /// (in serial)
    static std::map<std::size_t, std::vector<std::size_t>>
      contact_surface_map_volume_sweep_3d(Mesh& mesh, Function& u,
                                          const std::vector<std::size_t>& master_facets,
                                          const std::vector<std::size_t>& slave_facets);

  private:

    // Check whether two sets of tetrahedra collide
    static bool check_tet_set_collision(const Mesh& mmesh, std::size_t mi,
                                 const Mesh& smesh, std::size_t si);

    // Project surface forward from a facet using 'u', creating a prismoidal volume
    static std::vector<Point> create_deformed_segment_volume_3d(Mesh& mesh, std::size_t facet_index, const Function& u);

  };

}

#endif
