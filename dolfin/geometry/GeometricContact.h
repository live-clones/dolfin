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
    GeometricContact()
    {
      // Constructor
    }

    ~GeometricContact()
    {
      // Destructor
    }

    /// Calculate map from master facets to possible colliding slave facets
    void
      contact_surface_map_volume_sweep(Mesh& mesh, Function& u,
                                          const std::vector<std::size_t>& master_facets,
                                          const std::vector<std::size_t>& slave_facets);

    /// Get mapping
    const std::map<std::size_t, std::vector<std::size_t>>& master_to_slave() const
    {
      return _master_to_slave;
    }

  private:

    // Check whether two sets of triangles collide in 3D
    static bool check_tri_set_collision(const Mesh& mmesh, std::size_t mi,
                                        const Mesh& smesh, std::size_t si);

    // Check whether two sets of edges collide in 2D
    static bool check_edge_set_collision(const Mesh& mmesh, std::size_t mi,
                                        const Mesh& smesh, std::size_t si);

    // Project surface forward from a facet using 'u', creating a prismoidal volume in 2D or 3D
    static std::vector<Point> create_deformed_segment_volume(Mesh& mesh, std::size_t facet_index, const Function& u,
                                                             std::size_t gdim);

    std::map<std::size_t, std::vector<std::size_t>> _master_to_slave;

  };

}

#endif
