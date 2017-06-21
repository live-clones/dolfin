// Copyright (C) 2006 Anders Logg
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
// Modified by Garth N. Wells, 2008.
//
// First added:  2006-05-08
// Last changed: 2016-06-10

#ifndef __MESH_GEOMETRY_H
#define __MESH_GEOMETRY_H

#include <string>
#include <vector>
#include <dolfin/geometry/Point.h>
#include <dolfin/log/log.h>

namespace ufc
{
  class coordinate_mapping;
}

namespace dolfin
{
  class Function;

  /// MeshGeometry stores the geometry imposed on a mesh.

  /// Currently, the geometry is represented by the set of coordinates
  /// for the vertices of a mesh, but other representations are
  /// possible.
  class MeshGeometry
  {
  public:

    /// Create empty set of coordinates
    MeshGeometry();

    /// Copy constructor
    MeshGeometry(const MeshGeometry& geometry);

    /// Destructor
    ~MeshGeometry();

    /// Assignment
    const MeshGeometry& operator= (const MeshGeometry& geometry);

    /// Return Euclidean dimension of coordinate system
    std::size_t dim() const
    { return _dim; }

    /// Return polynomial degree of coordinate field
    std::size_t degree() const
    { return _degree; }

    /// Return the number of vertex coordinates
    std::size_t num_vertices() const
    {
      dolfin_assert(_coordinates.size() % _dim == 0);
      if (_degree > 1)
      {
        dolfin_assert(_entity_offsets.size() > 1);
        dolfin_assert(_entity_offsets[1].size() > 0);
        return _entity_offsets[1][0];
      }
      return _coordinates.size()/_dim;
    }

    /// Return the total number of points in the geometry, located on
    /// any entity
    std::size_t num_points() const
    {
      dolfin_assert(_coordinates.size() % _dim == 0);
      return _coordinates.size()/_dim;
    }

    /// Get vertex coordinates
    const double* vertex_coordinates(std::size_t point_index)
    {
      dolfin_assert(point_index < num_vertices());
      return &_coordinates[point_index*_dim];
    }

    /// Get vertex coordinates
    const double* point_coordinates(std::size_t point_index)
    {
      dolfin_assert(point_index*_dim < _coordinates.size());
      return &_coordinates[point_index*_dim];
    }

    /// Return value of coordinate with local index n in direction i
    double x(std::size_t n, std::size_t i) const
    {
      dolfin_assert((n*_dim + i) < _coordinates.size());
      dolfin_assert(i < _dim);
      return _coordinates[n*_dim + i];
    }

    /// Return array of values for coordinate with local index n
    const double* x(std::size_t n) const
    {
      dolfin_assert(n*_dim < _coordinates.size());
      return &_coordinates[n*_dim];
    }

    /// Return array of values for all coordinates
    std::vector<double>& x()
    { return _coordinates; }

    /// Return array of values for all coordinates
    const std::vector<double>& x() const
    { return _coordinates; }

    /// Return coordinate with local index n as a 3D point value
    Point point(std::size_t n) const;

    /// Initialize coordinate list to given dimension and degree
    void init(std::size_t dim, std::size_t degree);

    /// Initialise entities. To be called after init
    void init_entities(const std::vector<std::size_t>& num_entities);

    /// Get the number of coordinate points per entity for this degree
    std::size_t num_entity_coordinates(std::size_t entity_dim) const
    {
      // Calculate the number of points per entity for Lagrange
      // elements
      switch(entity_dim)
      {
      case 0:
        return 1;
      case 1:
        return (_degree - 1);
      case 2:
        return (_degree - 2)*(_degree - 1)/2;
      case 3:
        return (_degree - 3)*(_degree - 2)*(_degree - 1)/6;
      }
      dolfin_error("MeshGeometry.h",
                   "calculate number of points",
                   "Entity dimension out of range");
      return 0;
    }

    /// Get the index for an entity point in coordinates
    std::size_t get_entity_index(std::size_t entity_dim, std::size_t order,
                                 std::size_t index) const
    {
      dolfin_assert(entity_dim < _entity_offsets.size());
      dolfin_assert(order < _entity_offsets[entity_dim].size());
      const std::size_t idx = (_entity_offsets[entity_dim][order] + index);
      dolfin_assert(idx*_dim < _coordinates.size());
      return idx;
    }

    /// Set value of coordinate
    void set(std::size_t local_index, const double* x);

    /// Hash of coordinate values
    ///
    /// *Returns*
    ///     std::size_t
    ///         A tree-hashed value of the coordinates over all MPI processes
    ///
    std::size_t hash() const;

    /// Return informal string representation (pretty-print)
    std::string str(bool verbose) const;

    // NOTE: experimental and likely to change
    /// Set coordinate mapping for mesh
    void
      set_coordinate_mapping(std::shared_ptr<const ufc::coordinate_mapping> map)
    {
      _coordinate_mapping = map;
    }

    // NOTE: experimental and likely to change
    /// Return coordinate mapping (may be null)
    std::shared_ptr<const ufc::coordinate_mapping>
      get_coordinate_mapping() const
    {
      return _coordinate_mapping;
    }

  private:

    // Euclidean dimension
    std::size_t _dim;

    // Polynomial degree (1 = linear, 2 = quadratic etc.)
    std::size_t _degree;

    // Offsets to storage for coordinate points for each entity type
    std::vector<std::vector<std::size_t>> _entity_offsets;

    // Coordinates for all points stored as a contiguous array
    std::vector<double> _coordinates;

    // Coordinate mapping
    std::shared_ptr<const ufc::coordinate_mapping> _coordinate_mapping;

  };

}

#endif
