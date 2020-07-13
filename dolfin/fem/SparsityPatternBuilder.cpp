// Copyright (C) 2007-2010 Garth N. Wells
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
// Modified by Ola Skavhaug 2007
// Modified by Anders Logg 2008-2014
// Modified by Cecile Daversin-Catty 2018

#include <algorithm>

#include <dolfin/common/ArrayView.h>
#include <dolfin/common/MPI.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/log/log.h>
#include <dolfin/log/Progress.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MultiMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunctionSpace.h>
#include "MultiMeshDofMap.h"
#include "MultiMeshForm.h"
#include "SparsityPatternBuilder.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void
SparsityPatternBuilder::build(SparsityPattern& sparsity_pattern,
			      const Mesh& mesh,
			      const std::vector<const GenericDofMap*> dofmaps,
			      bool cells,
			      bool interior_facets,
			      bool exterior_facets,
			      bool vertices,
			      bool diagonal,
			      bool init,
			      bool finalize)
{
  /// Mono-domain version :
  /// Call the more generic build_mixed with mesh_ids[0] = mesh_ids[1] = integration mesh
  /// mesh_ids[0] is the mesh id associated with the test function
  /// mesh_ids[1] is the mesh id associated with the trial function
  std::vector<std::size_t> mesh_ids = {mesh.id(), mesh.id()};
  build_mixed(sparsity_pattern, mesh, mesh_ids, dofmaps, cells, interior_facets, exterior_facets, vertices, diagonal, init, finalize);
}
//-----------------------------------------------------------------------------
void
SparsityPatternBuilder::build_mixed(SparsityPattern& sparsity_pattern,
				    const Mesh& mesh,
				    std::vector<std::size_t> mesh_ids,
				    const std::vector<const GenericDofMap*> dofmaps,
				    bool cells,
				    bool interior_facets,
				    bool exterior_facets,
				    bool vertices,
				    bool diagonal,
				    bool init,
				    bool finalize)
{
  // Get global dimensions and local range
  const std::size_t rank = dofmaps.size();
  std::vector<std::shared_ptr<const IndexMap>> index_maps(rank);

  for (std::size_t i = 0; i < rank; ++i)
  {
    dolfin_assert(dofmaps[i]);
    index_maps[i] = dofmaps[i]->index_map();
  }

  // Initialise sparsity pattern
  if (init)
    sparsity_pattern.init(index_maps);

  // Only build for rank >= 2 (matrices and higher order tensors) that
  // require sparsity details
  if (rank < 2)
    return;

  // Vector to store macro-dofs, if required (for interior facets)
  std::vector<std::vector<dolfin::la_index>> macro_dofs(rank);

  // Create vector to point to dofs
  std::vector<ArrayView<const dolfin::la_index>> dofs(rank);
  std::vector<std::vector<dolfin::la_index>> dmaps(rank);

  // Build sparsity pattern for reals (globally supported basis members)
  // NOTE: It is very important that this is done before other integrals
  //       so that insertion of global nodes is no-op below
  // NOTE: We assume that global dofs contribute a whole row which is
  //       memory suboptimal (for restricted Lagrange multipliers) but very
  //       fast and certainly much better than quadratic scaling of usual
  //       insertion below
  std::vector<std::size_t> global_dofs0;
  dofmaps[sparsity_pattern.primary_dim()]->tabulate_global_dofs(global_dofs0);
  sparsity_pattern.insert_full_rows_local(global_dofs0);

  // FIXME: We iterate over the entire mesh even if the function space
  // is restricted. This works out fine since the local dofmap
  // returned on each cell will be an empty vector, but we might think
  // about optimizing this further.

  // Build sparsity pattern for cell integrals
  if (cells)
  {
    Progress p("Building sparsity pattern over cells", mesh.num_cells());
    auto mapping_map = mesh.topology().mapping();

    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      std::vector<std::vector<std::size_t>> cell_index(rank);
      std::vector<std::size_t> codim(rank);
      for (std::size_t i = 0; i < rank; ++i)
      {
	cell_index[i].push_back(cell->index());

	if(mesh_ids[i] != mesh.id() && mapping_map[mesh_ids[i]])
	{
	  auto mapping = mapping_map[mesh_ids[i]];
	  dolfin_assert(mapping->mesh()->id() == mesh_ids[i]);

	  codim[i] = mapping->mesh()->topology().dim() - mesh.topology().dim();
	  if(codim[i] == 0)
	    cell_index[i][0] = mapping->cell_map()[cell->index()];
	  else if(codim[i] == 1)
	  {
	    const std::size_t D = mapping->mesh()->topology().dim();
	    mapping->mesh()->init(D);
	    mapping->mesh()->init(D - 1, D);

	    Facet mesh_facet(*(mapping->mesh()), mapping->cell_map()[cell->index()]);
	    for(std::size_t j=0; j<mesh_facet.num_entities(D);j++)
	    {
	      Cell mesh_cell(*(mapping->mesh()), mesh_facet.entities(D)[j]);
	      if(j==0)
		cell_index[i][0] = mesh_cell.index();
	      else
		cell_index[i].push_back(mesh_cell.index());
	    }
	  }
#if 1 // 3D-1D work in progress #CHECKME
	  else if(codim[i] == 2)
          {
	    std::cout << "[SparsityBuilder] codim 2 - work in progress" << std::endl;
            const std::size_t D = mapping->mesh()->topology().dim();
            mapping->mesh()->init(D);
            mapping->mesh()->init(D - 2, D);

            Edge mesh_edge(*(mapping->mesh()), mapping->cell_map()[cell->index()]);
            for(std::size_t j=0; j<mesh_edge.num_entities(D);j++)
            {
              std::vector<double> cell_coordinates;
              Cell mesh_cell(*(mapping->mesh()), mesh_edge.entities(D)[j]);
              if(j==0)
                cell_index[i][0] = mesh_cell.index();
              else
                cell_index[i].push_back(mesh_cell.index());
            }
          }
#endif
	}
      }

      std::size_t nlocal_facets = cell_index[0].size();
      if(rank > 1)
	nlocal_facets = std::max(cell_index[0].size(), cell_index[1].size());

      for(std::size_t j=0; j<nlocal_facets; ++j)
      {
	for(std::size_t i=0; i<rank; ++i)
	{
	  std::size_t jidx  = (codim[i] != 0) ? j:0;
	  auto dmap = dofmaps[i]->cell_dofs(cell_index[i][jidx]);
	  dofs[i].set(dmap.size(), dmap.data());
	}
	sparsity_pattern.insert_local(dofs);
      }
      p++;
    }
  }

  // Build sparsity pattern for vertex/point integrals
  const std::size_t D = mesh.topology().dim();
  if (vertices)
  {
    mesh.init(0);
    mesh.init(0, D);

    std::vector<std::vector<dolfin::la_index>> global_dofs(rank);
    std::vector<std::vector<std::size_t>> local_to_local_dofs(rank);

    // Resize local dof map vector
    for (std::size_t i = 0; i < rank; ++i)
    {
      global_dofs[i].resize(dofmaps[i]->num_entity_dofs(0));
      local_to_local_dofs[i].resize(dofmaps[i]->num_entity_dofs(0));
    }

    Progress p("Building sparsity pattern over vertices", mesh.num_vertices());
    for (VertexIterator vert(mesh); !vert.end(); ++vert)
    {
      // Get mesh cell to which mesh vertex belongs (pick first)
      Cell mesh_cell(mesh, vert->entities(D)[0]);

      // Check that cell is not a ghost
      dolfin_assert(!mesh_cell.is_ghost());

      // Get local index of vertex with respect to the cell
      const std::size_t local_vertex = mesh_cell.index(*vert);
      for (std::size_t i = 0; i < rank; ++i)
      {
        auto dmap = dofmaps[i]->cell_dofs(mesh_cell.index());
        dofs[i].set(dmap.size(), dmap.data());
        dofmaps[i]->tabulate_entity_dofs(local_to_local_dofs[i], 0,
                                         local_vertex);

        // Copy cell dofs to local dofs and tabulated values to
        for (std::size_t j = 0; j < local_to_local_dofs[i].size(); ++j)
          global_dofs[i][j] = dofs[i][local_to_local_dofs[i][j]];
      }

      // Insert non-zeroes in sparsity pattern
      std::vector<ArrayView<const dolfin::la_index>> global_dofs_p(rank);
      for (std::size_t i = 0; i < rank; ++i)
        global_dofs_p[i].set(global_dofs[i]);
      sparsity_pattern.insert_local(global_dofs_p);
      p++;
    }
  }

  // Note: no need to iterate over exterior facets since those dofs
  //       are included when tabulating dofs on all cells

  // Build sparsity pattern for interior/exterior facet integrals
  if (interior_facets || exterior_facets)
  {
    // Compute facets and facet - cell connectivity if not already
    // computed
    mesh.init(D - 1);
    mesh.init(D - 1, D);
    if (!mesh.ordered())
    {
      dolfin_error("SparsityPatternBuilder.cpp",
                   "compute sparsity pattern",
                   "Mesh is not ordered according to the UFC numbering convention. "
                   "Consider calling mesh.order()");
    }

    Progress p("Building sparsity pattern over interior facets",
               mesh.num_facets());
    for (FacetIterator facet(mesh); !facet.end(); ++facet)
    {
      bool this_exterior_facet = false;
      if (facet->num_global_entities(D) == 1)
        this_exterior_facet = true;

      // Check facet type
      if (exterior_facets && this_exterior_facet && !cells)
      {
        // Get cells incident with facet
        dolfin_assert(facet->num_entities(D) == 1);
        Cell cell(mesh, facet->entities(D)[0]);

        // Tabulate dofs for each dimension and get local dimensions
        for (std::size_t i = 0; i < rank; ++i)
        {
          auto dmap = dofmaps[i]->cell_dofs(cell.index());
          dofs[i].set(dmap.size(), dmap.data());
        }

        // Insert dofs
        sparsity_pattern.insert_local(dofs);
      }
      else if (interior_facets && !this_exterior_facet)
      {
        if (facet->num_entities(D) == 1)
        {
          dolfin_assert(facet->is_ghost());
          continue;
        }

        // Get cells incident with facet
        dolfin_assert(facet->num_entities(D) == 2);
        Cell cell0(mesh, facet->entities(D)[0]);
        Cell cell1(mesh, facet->entities(D)[1]);

        // Tabulate dofs for each dimension on macro element
        for (std::size_t i = 0; i < rank; i++)
        {
          // Get dofs for each cell
          auto cell_dofs0 = dofmaps[i]->cell_dofs(cell0.index());
          auto cell_dofs1 = dofmaps[i]->cell_dofs(cell1.index());

          // Create space in macro dof vector
          macro_dofs[i].resize(cell_dofs0.size() + cell_dofs1.size());

          // Copy cell dofs into macro dof vector
          std::copy(cell_dofs0.data(), cell_dofs0.data() + cell_dofs0.size(),
                    macro_dofs[i].begin());
          std::copy(cell_dofs1.data(), cell_dofs1.data() + cell_dofs1.size(),
                    macro_dofs[i].begin() + cell_dofs0.size());

          // Store pointer to macro dofs
          dofs[i].set(macro_dofs[i]);
        }

        // Insert dofs
        sparsity_pattern.insert_local(dofs);
      }
      p++;
    }
  }

  if (diagonal)
  {
    dolfin_assert(rank == 2);
    const std::size_t primary_dim = sparsity_pattern.primary_dim();
    const std::size_t primary_codim = primary_dim == 0 ? 1 : 0;
    const std::pair<std::size_t, std::size_t> primary_range
      = index_maps[primary_dim]->local_range();
    const std::size_t secondary_range
      = index_maps[primary_codim]->size(IndexMap::MapSize::GLOBAL);
    const std::size_t diagonal_range
      = std::min(primary_range.second, secondary_range);

    Progress p("Building sparsity pattern over diagonal",
               diagonal_range - primary_range.first);
    for (std::size_t j = primary_range.first; j < diagonal_range; j++)
    {
      sparsity_pattern.insert_global(j, j);
      p++;
    }
  }

  // Finalize sparsity pattern (communicate off-process terms)
  if (finalize)
    sparsity_pattern.apply();
}
//-----------------------------------------------------------------------------
void SparsityPatternBuilder::build_multimesh_sparsity_pattern(
  SparsityPattern& sparsity_pattern,
  const MultiMeshForm& form)
{
  // Get global dimensions and local range
  const std::size_t rank = form.rank();
  std::vector<std::shared_ptr<const IndexMap>> index_maps(rank);
  for (std::size_t i = 0; i < rank; ++i)
    index_maps[i] = form.function_space(i)->dofmap()->index_map();

  // Initialize sparsity pattern
  sparsity_pattern.init(index_maps);

  // Iterate over each part
  for (std::size_t part = 0; part < form.num_parts(); part++)
  {
    // Get mesh on current part (assume it's the same for all
    // arguments)
    const Mesh& mesh = *form.function_space(0)->part(part)->mesh();

    // Build list of dofmaps
    std::vector<const GenericDofMap*> dofmaps;
    for (std::size_t i = 0; i < form.rank(); i++)
      dofmaps.push_back(&*form.function_space(i)->dofmap()->part(part));

    log(PROGRESS, "Building intra-mesh sparsity pattern on part %d.", part);

    // Build sparsity pattern for part by calling the regular dofmap
    // builder. This builds the sparsity pattern for all interacting
    // dofs on the current part.
    build(sparsity_pattern, mesh, dofmaps,
	  true, false, false, true, false, false);

    log(PROGRESS, "Building inter-mesh sparsity pattern on part %d.", part);

    // Build sparsity pattern for interface. This builds the sparsity
    // pattern for all dofs that may interact across the interface
    // between cutting meshes.
    _build_multimesh_sparsity_pattern_interface(sparsity_pattern, form, part);
  }

  log(PROGRESS, "Applying changes to sparsity pattern.");

  // Finalize sparsity pattern
  sparsity_pattern.apply();
}
//-----------------------------------------------------------------------------
void SparsityPatternBuilder::_build_multimesh_sparsity_pattern_interface(
  SparsityPattern& sparsity_pattern,
  const MultiMeshForm& form,
  std::size_t part)
{
  // Get multimesh
  const auto& multimesh = form.multimesh();

  // Get collision map
  const auto& cmap = multimesh->collision_map_cut_cells(part);

  // Data structures for storing dofs on cut (0) and cutting cell (1)
  std::vector<ArrayView<const dolfin::la_index>> dofs_0(form.rank());
  std::vector<ArrayView<const dolfin::la_index>> dofs_1(form.rank());

  // FIXME: We need two different lists here because the interface
  // FIXME: of insert() requires a list of pointers to dofs. Consider
  // FIXME: improving the interface of SparsityPattern.

  // Data structure for storing dofs on macro cell (0 + 1)
  std::vector<std::vector<dolfin::la_index>> dofs(form.rank());
  std::vector<ArrayView<const dolfin::la_index>> _dofs(form.rank());

  // Iterate over all cut cells in collision map
  for (auto it = cmap.begin(); it != cmap.end(); ++it)
  {
    // Get cut cell index
    const unsigned int cut_cell_index = it->first;

    // Get dofs for cut cell
    for (std::size_t i = 0; i < form.rank(); i++)
    {
      const auto& dofmap = form.function_space(i)->dofmap()->part(part);
      auto dmap = dofmap->cell_dofs(cut_cell_index);
      dofs_0[i].set(dmap.size(), dmap.data());
    }

    // Iterate over cutting cells
    const auto& cutting_cells = it->second;
    for (auto jt = cutting_cells.begin(); jt != cutting_cells.end(); jt++)
    {
      // Get cutting part and cutting cell index
      const std::size_t cutting_part = jt->first;
      const std::size_t cutting_cell_index = jt->second;

      // Add dofs for cutting cell
      for (std::size_t i = 0; i < form.rank(); i++)
      {
        // Get dofs for cutting cell
        const auto& dofmap
          = form.function_space(i)->dofmap()->part(cutting_part);
        auto dmap = dofmap->cell_dofs(cutting_cell_index);
        dofs_1[i].set(dmap.size(), dmap.data());

        // Collect dofs for cut and cutting cell
        dofs[i].resize(dofs_0[i].size() + dofs_1[i].size());
        std::copy(dofs_0[i].begin(), dofs_0[i].end(), dofs[i].begin());
        std::copy(dofs_1[i].begin(), dofs_1[i].end(),
                  dofs[i].begin() + dofs_0[i].size());
        _dofs[i].set(dofs[i]);
      }

      // Insert into sparsity pattern
      sparsity_pattern.insert_local(_dofs);
    }
  }
}
//-----------------------------------------------------------------------------
