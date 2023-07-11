// Copyright (C) 2006-2016 Anders Logg
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
// Modified by Johan Hoffman 2007
// Modified by Garth N. Wells 2007-2011
// Modified by Niclas Jansson 2008
// Modified by Kristoffer Selim 2008
// Modified by Andre Massing 2009-2010
// Modified by Johannes Ring 2012
// Modified by Marie E. Rognes 2012
// Modified by Mikael Mortensen 2012
// Modified by Jan Blechta 2013
// Modified by Cecile Daversin-Catty 2018
//
// First added:  2006-05-09
// Last changed: 2019-12-04

#include <dolfin/ale/ALE.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/MPI.h>
#include <dolfin/common/Timer.h>
#include <dolfin/common/utils.h>
#include <dolfin/function/Expression.h>
#include <dolfin/io/File.h>
#include <dolfin/log/log.h>
#include <dolfin/geometry/BoundingBoxTree.h>
#include "BoundaryMesh.h"
#include "Cell.h"
#include "DistributedMeshTools.h"
#include "Facet.h"
#include "LocalMeshData.h"
#include "MeshColoring.h"
#include "MeshOrdering.h"
#include "MeshPartitioning.h"
#include "MeshRenumbering.h"
#include "MeshSmoothing.h"
#include "MeshTransformation.h"
#include "TopologyComputation.h"
#include "Vertex.h"
#include "Mesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Mesh::Mesh() : Mesh(MPI_COMM_WORLD)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(MPI_Comm comm) : Variable("mesh", "DOLFIN mesh"),
                            Hierarchical<Mesh>(*this), _ordered(false),
                            _mpi_comm(comm), _ghost_mode("none")
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(const Mesh& mesh) : Variable("mesh", "DOLFIN mesh"),
                               Hierarchical<Mesh>(*this), _ordered(false),
                               _mpi_comm(mesh.mpi_comm()),
                               _ghost_mode("none")
{
  *this = mesh;
}
//-----------------------------------------------------------------------------
Mesh::Mesh(std::string filename) : Mesh(MPI_COMM_WORLD, filename)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Mesh::Mesh(MPI_Comm comm, std::string filename)
  : Variable("mesh", "DOLFIN mesh"), Hierarchical<Mesh>(*this), _ordered(false),
  _mpi_comm(comm), _ghost_mode("none")
{
  File file(_mpi_comm.comm(), filename);
  file >> *this;
}
//-----------------------------------------------------------------------------
Mesh::Mesh(MPI_Comm comm, LocalMeshData& local_mesh_data)
  : Variable("mesh", "DOLFIN mesh"), Hierarchical<Mesh>(*this),
  _ordered(false), _mpi_comm(comm), _ghost_mode("none")
{
  const std::string ghost_mode = parameters["ghost_mode"];
  MeshPartitioning::build_distributed_mesh(*this, local_mesh_data, ghost_mode);
}
//-----------------------------------------------------------------------------
Mesh::~Mesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
const Mesh& Mesh::operator=(const Mesh& mesh)
{
  // Assign data
  _topology = mesh._topology;
  _geometry = mesh._geometry;
  _domains = mesh._domains;
  _data = mesh._data;
  if (mesh._cell_type)
    _cell_type.reset(CellType::create(mesh._cell_type->cell_type()));
  else
    _cell_type.reset();
  _ordered = mesh._ordered;
  _cell_orientations = mesh._cell_orientations;
  _ghost_mode = mesh._ghost_mode;

  // Rename
  rename(mesh.name(), mesh.label());

  // Call assignment operator for base class
  Hierarchical<Mesh>::operator=(mesh);

  return *this;
}
//-----------------------------------------------------------------------------
MeshData& Mesh::data()
{
  return _data;
}
//-----------------------------------------------------------------------------
const MeshData& Mesh::data() const
{
  return _data;
}
//-----------------------------------------------------------------------------
std::size_t Mesh::init(std::size_t dim) const
{
  // This function is obviously not const since it may potentially
  // compute new connectivity. However, in a sense all connectivity of
  // a mesh always exists, it just hasn't been computed yet. The
  // const_cast is also needed to allow iterators over a const Mesh to
  // create new connectivity.

  // Skip if mesh is empty
  if (num_cells() == 0)
  {
    //Do not display empty mesh warning if we are dealing with submeshes
    if(_topology.mapping().empty())
      warning("Mesh is empty, unable to create entities of dimension %d.", dim);
    return 0;
  }

  // Skip if already computed
  if (_topology.size(dim) > 0)
    return _topology.size(dim);

  // Skip vertices and cells (should always exist)
  if (dim == 0 || dim == _topology.dim())
    return _topology.size(dim);

  // Check that mesh is ordered
  if (!ordered())
  {
    dolfin_error("Mesh.cpp",
                 "initialize mesh entities",
                 "Mesh is not ordered according to the UFC numbering convention. Consider calling mesh.order()");
  }

  // Compute connectivity
  Mesh* mesh = const_cast<Mesh*>(this);
  TopologyComputation::compute_entities(*mesh, dim);

  // Order mesh if necessary
  if (!ordered())
    mesh->order();

  return _topology.size(dim);
}
//-----------------------------------------------------------------------------
void Mesh::init(std::size_t d0, std::size_t d1) const
{
  // This function is obviously not const since it may potentially
  // compute new connectivity. However, in a sense all connectivity of
  // a mesh always exists, it just hasn't been computed yet. The
  // const_cast is also needed to allow iterators over a const Mesh to
  // create new connectivity.

  // Skip if mesh is empty
  if (num_cells() == 0)
  {
    //Do not display empty mesh warning if we are dealing with submeshes
    if(_topology.mapping().empty())
      warning("Mesh is empty, unable to create connectivity %d --> %d.", d0, d1);
    return;
  }

  // Skip if already computed
  if (!_topology(d0, d1).empty())
    return;

  // Check that mesh is ordered
  if (!ordered())
  {
    dolfin_error("Mesh.cpp",
                 "initialize mesh connectivity",
                 "Mesh is not ordered according to the UFC numbering convention. Consider calling mesh.order()");
  }

  // Compute connectivity
  Mesh* mesh = const_cast<Mesh*>(this);
  TopologyComputation::compute_connectivity(*mesh, d0, d1);

  // Order mesh if necessary
  if (!ordered())
    mesh->order();
}
//-----------------------------------------------------------------------------
void Mesh::init() const
{
  // Compute all entities
  for (std::size_t d = 0; d <= topology().dim(); d++)
    init(d);

  // Compute all connectivity
  for (std::size_t d0 = 0; d0 <= topology().dim(); d0++)
    for (std::size_t d1 = 0; d1 <= topology().dim(); d1++)
      init(d0, d1);
}
//-----------------------------------------------------------------------------
void Mesh::init_global(std::size_t dim) const
{
  init(dim);
  DistributedMeshTools::number_entities(*this, dim);
}
//-----------------------------------------------------------------------------
void Mesh::clean()
{
  const std::size_t D = topology().dim();
  for (std::size_t d0 = 0; d0 <= D; d0++)
  {
    for (std::size_t d1 = 0; d1 <= D; d1++)
    {
      if (!(d0 == D && d1 == 0))
        _topology.clear(d0, d1);
    }
  }
}
//-----------------------------------------------------------------------------
void Mesh::order()
{
  // Order mesh
  MeshOrdering::order(*this);

  // Remember that the mesh has been ordered
  _ordered = true;

  // Clear any cell_orientations (as these depend on the ordering)
  _cell_orientations.clear();
}
//-----------------------------------------------------------------------------
bool Mesh::ordered() const
{
  // Don't check if we know (or think we know) that the mesh is ordered
  if (_ordered)
    return true;

  _ordered = MeshOrdering::ordered(*this);
  return _ordered;
}
//-----------------------------------------------------------------------------
dolfin::Mesh Mesh::renumber_by_color() const
{
  const std::size_t D = topology().dim();
  const std::vector<std::size_t> coloring_type = {{D, 0, D}};
  return MeshRenumbering::renumber_by_color(*this, coloring_type);
}
//-----------------------------------------------------------------------------
void Mesh::scale(double factor)
{
  MeshTransformation::scale(*this, factor);
}
//-----------------------------------------------------------------------------
void Mesh::translate(const Point& point)
{
  MeshTransformation::translate(*this, point);
}
//-----------------------------------------------------------------------------
void Mesh::rotate(double angle, std::size_t axis)
{
  MeshTransformation::rotate(*this, angle, axis);
}
//-----------------------------------------------------------------------------
void Mesh::rotate(double angle, std::size_t axis, const Point& point)
{
  MeshTransformation::rotate(*this, angle, axis, point);
}
//-----------------------------------------------------------------------------
void Mesh::smooth(std::size_t num_iterations)
{
  MeshSmoothing::smooth(*this, num_iterations);
}
//-----------------------------------------------------------------------------
void Mesh::smooth_boundary(std::size_t num_iterations, bool harmonic_smoothing)
{
  MeshSmoothing::smooth_boundary(*this, num_iterations, harmonic_smoothing);
}
//-----------------------------------------------------------------------------
void Mesh::snap_boundary(const SubDomain& sub_domain, bool harmonic_smoothing)
{
  MeshSmoothing::snap_boundary(*this, sub_domain, harmonic_smoothing);
}
//-----------------------------------------------------------------------------
const std::vector<std::size_t>& Mesh::color(std::string coloring_type) const
{
  // Define graph type
  const std::size_t dim = MeshColoring::type_to_dim(coloring_type, *this);
  const std::vector<std::size_t> _coloring_type
    = {{topology().dim(), dim, topology().dim()}};

  return color(_coloring_type);
}
//-----------------------------------------------------------------------------
const std::vector<std::size_t>&
Mesh::color(std::vector<std::size_t> coloring_type) const
{
  // Find color data
  std::map<std::vector<std::size_t>, std::pair<std::vector<std::size_t>,
           std::vector<std::vector<std::size_t>>>>::const_iterator
    coloring_data = this->topology().coloring.find(coloring_type);

  if (coloring_data != this->topology().coloring.end())
  {
    dolfin_debug("Mesh has already been colored, not coloring again.");
    return coloring_data->second.first;
  }

  // We do the same const-cast trick here as in the init() functions
  // since we are not really changing the mesh, just attaching some
  // auxiliary data to it.
  Mesh* _mesh = const_cast<Mesh*>(this);
  return MeshColoring::color(*_mesh, coloring_type);
}
//-----------------------------------------------------------------------------
std::shared_ptr<BoundingBoxTree> Mesh::bounding_box_tree() const
{
  // Allocate and build tree if necessary
  if (!_tree)
  {
    _tree.reset(new BoundingBoxTree());
    _tree->build(*this);
  }

  return _tree;
}
//-----------------------------------------------------------------------------
double Mesh::hmin() const
{
  double h = std::numeric_limits<double>::max();
  for (CellIterator cell(*this); !cell.end(); ++cell)
    h = std::min(h, cell->h());

  return h;
}
//-----------------------------------------------------------------------------
double Mesh::hmax() const
{
  double h = 0.0;
  for (CellIterator cell(*this); !cell.end(); ++cell)
    h = std::max(h, cell->h());

  return h;
}
//-----------------------------------------------------------------------------
double Mesh::rmin() const
{
  double r = std::numeric_limits<double>::max();
  for (CellIterator cell(*this); !cell.end(); ++cell)
    r = std::min(r, cell->inradius());

  return r;
}
//-----------------------------------------------------------------------------
double Mesh::rmax() const
{
  double r = 0.0;
  for (CellIterator cell(*this); !cell.end(); ++cell)
    r = std::max(r, cell->inradius());

  return r;
}
//-----------------------------------------------------------------------------
std::size_t Mesh::hash() const
{
  // Get local hashes
  const std::size_t kt_local = _topology.hash();
  const std::size_t kg_local = _geometry.hash();

  // Compute global hash
  const std::size_t kt = hash_global(_mpi_comm.comm(), kt_local);
  const std::size_t kg = hash_global(_mpi_comm.comm(), kg_local);

  // Compute hash based on the Cantor pairing function
  return (kt + kg)*(kt + kg + 1)/2 + kg;
}
//-----------------------------------------------------------------------------
std::string Mesh::str(bool verbose) const
{
  std::stringstream s;
  if (verbose)
  {
    s << str(false) << std::endl << std::endl;

    s << indent(_geometry.str(true));
    s << indent(_topology.str(true));
    s << indent(_data.str(true));
  }
  else
  {
    std::string cell_type("undefined cell type");
    if (_cell_type)
      cell_type = _cell_type->description(true);

    s << "<Mesh of topological dimension "
      << topology().dim() << " ("
      << cell_type << ") with "
      << num_vertices() << " vertices and "
      << num_cells() << " cells, "
      << (_ordered ? "ordered" : "unordered") << ">";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
const std::vector<int>& Mesh::cell_orientations() const
{
  return _cell_orientations;
}
//-----------------------------------------------------------------------------
void Mesh::init_cell_orientations(const Expression& global_normal)
{
  std::size_t gdim = geometry().dim();
  std::size_t ndim = global_normal.value_size();

  // Check that global_normal has the "right" size
  // Allowing 3 if gdim < 3 to avoid breaking legacy code.
  if (ndim < gdim && ndim <= 3)
  {
     dolfin_error("Mesh.cpp",
                  "initialize cell orientations",
                  "Global normal value size is %d, smaller than gdim (%d)",
                  ndim, gdim);
  }

  // Resize storage
  _cell_orientations.resize(num_cells());

  // Set orientation
  Array<double> values(ndim);
  Point up;
  for (CellIterator cell(*this); !cell.end(); ++cell)
  {
    // Extract cell midpoint as Array
    const Array<double> x(3, cell->midpoint().coordinates());

    // Evaluate global normal at cell midpoint
    global_normal.eval(values, x);

    // Extract values as Point
    for (unsigned int i = 0; i < ndim; i++)
      up[i] = values[i];
    for (unsigned int i = ndim; i < gdim; i++)
      up[i] = 0.0;

    // Set orientation as orientation relative to up direction.
    dolfin_assert(cell->index() < _cell_orientations.size());
    _cell_orientations[cell->index()] = cell->orientation(up);
  }
}
//-----------------------------------------------------------------------------
std::string Mesh::ghost_mode() const
{
  dolfin_assert(_ghost_mode == "none"
                || _ghost_mode == "shared_vertex"
                || _ghost_mode == "shared_facet");
  return _ghost_mode;
}
//-----------------------------------------------------------------------------

void Mesh::build_mapping(std::shared_ptr<const Mesh> other) const
{
  // Find a parent mesh shared by the two meshes
  std::vector<unsigned> current_keys, other_keys, common_keys;
  for(auto mapping : this->_topology.mapping())
    current_keys.push_back(mapping.first);
  for(auto mapping : other->_topology.mapping())
    other_keys.push_back(mapping.first);
  std::set_intersection(current_keys.begin(),current_keys.end(),other_keys.begin(),other_keys.end(),back_inserter(common_keys));

  if (common_keys.size() == 0)
    throw std::logic_error("Cannot find common parent mesh");
  unsigned parent_id = common_keys[0];

  // NOTE : This function is only used in MixedAssembler
  // for now, and needs only the cell_map to be filled.
  // We might need to fill the corresponding vertex_map
  // at some point.
  std::vector<std::size_t> current_cell_map = this->topology().mapping()[parent_id]->cell_map();
  std::vector<std::size_t> other_cell_map = other->topology().mapping()[parent_id]->cell_map();
  std::vector<std::size_t> new_vertex_map;
  std::vector<std::size_t> new_cell_map;
  new_cell_map.reserve(current_cell_map.size());

  // Check if child mesh has no cells on process
  bool empty_child_mesh = (other_cell_map.size() == 0);

  std::vector<std::size_t>::iterator found;
  // Build the new mapping current <-> other from the common parent (cell_map only)

  for (std::size_t i = 0; i < current_cell_map.size(); ++i)
  {
    bool new_idx_found = false;
    // Co-dimension 0
    if (this->topology().dim() == other->topology().dim())
    {
      found = std::find(other_cell_map.begin(), other_cell_map.end(), current_cell_map[i]);
      if (found != other_cell_map.end())
      {
        new_idx_found = true;
        new_cell_map.push_back(found - other_cell_map.begin());
      }
    }
    // Co-dimension 1
    else if (this->topology().dim() == other->topology().dim() - 1)
    {
      // current_cell_map[i] is the index of a facet in the parent mesh.
      // We need to find the entities of dimension dim(other) = dim(current) + 1
      // related to this facet.

      // mesh_facet is the facet of the parent mesh corresponding to current_cell_map[i]
      Facet mesh_facet(*(this->topology().mapping()[parent_id]->mesh()), current_cell_map[i]);
      auto D = other->topology().dim();
      this->topology().mapping()[parent_id]->mesh()->init(D);
      this->topology().mapping()[parent_id]->mesh()->init(D - 1, D);
      other->init(this->topology().dim());
      // Find a cell in <other> owning mesh_facet
      for (std::size_t j = 0; j < mesh_facet.num_entities(D); j++)
      {
        Cell mesh_cell(*(this->topology().mapping()[parent_id]->mesh()), mesh_facet.entities(other->topology().dim())[j]);
        found = std::find(other_cell_map.begin(), other_cell_map.end(), mesh_cell.index());
        if (found != other_cell_map.end())
        {
          auto position = found - other_cell_map.begin();
          Cell other_mesh_cell(*(other), position);
          // Find which facet of this cell is the one we are looking for
          for (std::size_t k = 0; k < other_mesh_cell.num_entities(this->topology().dim()); k++)
          {
            // NOTE : Comparison based on midpoint coordinates
            // Is it still ok with very small mesh size ?
            Facet other_mesh_facet(*(other), other_mesh_cell.entities(this->topology().dim())[k]);
            if (other_mesh_facet.midpoint() == mesh_facet.midpoint())
            {
              new_idx_found = true;
              new_cell_map.push_back(other_mesh_cell.entities(this->topology().dim())[k]);
              break;
            }
          }
          break;
        }
      }
    }
    else
     dolfin_error("Mesh.cpp",
                  "build_mapping",
                  "The dimension of the mesh given as parameter (%d) cannot be lower than the dimension of the current mesh (%d)",
                  other->topology().dim(), this->topology().dim());

    if ((!empty_child_mesh) && (!new_idx_found))
    {
      std::cout << "Error in building the mapping ("
                << this->id() << ", " << other->id() << ") :"
                << "Index not found." << std::endl;
    }
  }
  this->_topology.add_mapping(std::make_pair(other->id(), std::make_shared<MeshView>(other, new_vertex_map, new_cell_map)));
}
