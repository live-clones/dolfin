// Copyright (C) 2013-2016 Garth N. Wells, Steven Vandekerckhove,
// Tormod Landet, Jan Blechta
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

#include <array>
#include <memory>
#include <vector>
#include <Eigen/Dense>

#include <dolfin/common/ArrayView.h>
#include <dolfin/common/Set.h>
#include <dolfin/common/Timer.h>
#include <dolfin/common/types.h>
#include <dolfin/fem/LocalAssembler.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericLinearAlgebraFactory.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/log/log.h>
#include <dolfin/log/Progress.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include "assemble.h"
#include "Form.h"
#include "GenericDofMap.h"
#include "UFC.h"
#include "LocalPatchSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
LocalPatchSolver::LocalPatchSolver(std::vector<std::shared_ptr<const Form>> a,
                                   std::vector<std::shared_ptr<const Form>> L,
                                   SolverType solver_type,
                                   BCType bc_type)
  : _a(a), _formL(L), _solver_type(solver_type), _bc_type(bc_type),
    _num_patches(a.size())
{
  _check_input(&a, &L);
}
//-----------------------------------------------------------------------------
LocalPatchSolver::LocalPatchSolver(std::vector<std::shared_ptr<const Form>> a,
                                   SolverType solver_type,
                                   BCType bc_type)
  : _a(a), _solver_type(solver_type), _bc_type(bc_type),
    _num_patches(a.size())
{
  _check_input(&a, nullptr);
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::solve_global_rhs(Function& u) const
{
  solve_global_rhs(std::vector<Function*>(_num_patches, &u));
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::solve_global_rhs(std::vector<Function*> u) const
{
  _check_input(u);

  // Compute RHSs (global)
  std::vector<std::shared_ptr<GenericVector>> b;
  for (std::size_t i = 0; i < _num_patches; ++i)
  {
    b.push_back(u[i]->vector()->factory().create_vector(u[i]->vector()->mpi_comm()));
    dolfin_assert(_formL[i]);
    assemble(*b[i], *_formL[i]);
  }

  // Extract the vectors where the solution will be stored
  std::vector<GenericVector*> x;
  for (auto ui: u)
  {
    dolfin_assert(ui);
    x.push_back(ui->vector().get());
  }

  // Solve local problems
  std::vector<const GenericVector*> b_plain;
  for (const auto bi: b)
    b_plain.push_back(bi.get());
  _solve_local(x, &b_plain, nullptr);
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::solve_local_rhs(Function& u) const
{
  solve_local_rhs(std::vector<Function*>(_num_patches, &u));
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::solve_local_rhs(std::vector<Function*> u) const
{
  _check_input(u);

  // Extract the vectors where the solution will be stored
  std::vector<GenericVector*> x;
  for (auto ui: u)
  {
    dolfin_assert(ui);
    x.push_back(ui->vector().get());
  }

  // Loop over all cells and assemble local LHS & RHS which are then solved
  _solve_local(x, nullptr, nullptr);
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::solve_local(GenericVector& x, const GenericVector& b,
                                   const GenericDofMap& dofmap_b) const
{
  const std::vector<GenericVector*> _x(_num_patches, &x);
  const std::vector<const GenericVector*> _b(_num_patches, &b);
  std::vector<const GenericDofMap*> _dofmap_b(_num_patches, &dofmap_b);
  _solve_local(_x, &_b, &_dofmap_b);
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::solve_local(std::vector<GenericVector*> x,
                                   std::vector<const GenericVector*> b,
                                   std::vector<const GenericDofMap*> dofmap_b) const
{
  _check_input(x, b, dofmap_b);
  _solve_local(x, &b, &dofmap_b);
}
//----------------------------------------------------------------------------
// FIXME: Should the function be in a namespace? It's not in any header.
// FIXME: It is suboptimal to: first, merge cell dofs into patch dofs;
//        second, sort them; third, reconstruct a map using this function
inline std::vector<dolfin::la_index> _compute_cell_to_patch_map(
  const ArrayView<const dolfin::la_index> cell_dofs,
  const ArrayView<const dolfin::la_index> patch_dofs)
{
  std::vector<dolfin::la_index> map;
  map.reserve(cell_dofs.size());
  for (const dolfin::la_index dof: cell_dofs)
  {
    const auto it = std::find(patch_dofs.begin(), patch_dofs.end(), dof);
    if (it != patch_dofs.end())
      map.push_back(std::distance(patch_dofs.begin(), it));
  }
  return map;
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::_solve_local(std::vector<GenericVector*> x,
  const std::vector<const GenericVector*>* const global_b,
  std::vector<const GenericDofMap*>* dofmap_L) const
{
  // Set timer
  Timer timer("Solve local patch problems");

  // Create UFC objects
  std::vector<UFC> ufc_a;
  for (auto a: _a)
    ufc_a.emplace_back(*a);
  std::vector<UFC> ufc_L;

  std::vector<const GenericDofMap*> _dofmap_L_vec;

  // Check that we have valid linear form or a dofmap for it
  if (dofmap_L)
    dolfin_assert(global_b);
  else
  {
    dofmap_L = &_dofmap_L_vec;
    for (auto L: _formL)
    {
      dolfin_assert(L);
      dolfin_assert(L->rank() == 1);
      dolfin_assert(L->function_space(0)->dofmap());
      _dofmap_L_vec.push_back(L->function_space(0)->dofmap().get());
      ufc_L.emplace_back(*L);
    }
  }

  // Extract the mesh
  dolfin_assert(_a[0]->function_space(0)->mesh());
  const Mesh& mesh = *_a[0]->function_space(0)->mesh();

  // Get bilinear form dofmaps
  std::array<std::shared_ptr<const GenericDofMap>, 2> dofmaps_a
    = {{_a[0]->function_space(0)->dofmap(), _a[0]->function_space(1)->dofmap()}};
  dolfin_assert(dofmaps_a[0] && dofmaps_a[1]);

  // Extract domains of forms
  std::vector<const MeshFunction<std::size_t>*> cell_domains_a;
  std::vector<const MeshFunction<std::size_t>*> exterior_facet_domains_a;
  std::vector<const MeshFunction<std::size_t>*> interior_facet_domains_a;
  std::vector<const MeshFunction<std::size_t>*> cell_domains_L;
  std::vector<const MeshFunction<std::size_t>*> exterior_facet_domains_L;
  std::vector<const MeshFunction<std::size_t>*> interior_facet_domains_L;
  for (const auto a: _a)
  {
    cell_domains_a.push_back(a->cell_domains().get());
    exterior_facet_domains_a.push_back(a->exterior_facet_domains().get());
    interior_facet_domains_a.push_back(a->interior_facet_domains().get());
  }
  if (!global_b)
  {
    for (const auto L: _formL)
    {
      cell_domains_L.push_back(L->cell_domains().get());
      exterior_facet_domains_L.push_back(L->exterior_facet_domains().get());
      interior_facet_domains_L.push_back(L->interior_facet_domains().get());
    }
  }

  // Eigen data structures and factorisations for cell data structures
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    A_e, b_e, A_cell, b_cell;
  Eigen::VectorXd x_e;
  Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                    Eigen::RowMajor>> lu;
  Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                           Eigen::RowMajor>> cholesky;
  const bool use_cache = !(_cholesky_cache.empty() && _lu_cache.empty());

  // Zero tensors
  for (auto xi: x)
  {
    dolfin_assert(xi);
    xi->zero();
  }

  // Loop over vertices and solve local patch problems
  Progress p("Performing local (patch-wise) solve", mesh.num_vertices());
  Set<la_index> dofs_a0, dofs_a1, dofs_L;
  ArrayView<const la_index> dofs_a0_ptr, dofs_a1_ptr, dofs_L_ptr;
  // FIXME: Is the iterator correctly ghosted?
  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
  {
    // Clear the patch dofs
    dofs_a0.clear();
    dofs_a1.clear();
    dofs_L.clear();

    // Get local dofs cell-by-cell on the patch
    for (CellIterator cell(*vertex); !cell.end(); ++cell)
    {
      // Get local-to-global dof maps for cell
      const ArrayView<const dolfin::la_index> celldofs_a0
        = dofmaps_a[0]->cell_dofs(cell->index());
      const ArrayView<const dolfin::la_index> celldofs_a1
        = dofmaps_a[1]->cell_dofs(cell->index());
      const ArrayView<const dolfin::la_index> celldofs_L
        = (*dofmap_L)[0]->cell_dofs(cell->index());

      // Add cell dofs to patch dofs list
      dofs_a0.insert(celldofs_a0.begin(), celldofs_a0.end());
      dofs_a1.insert(celldofs_a1.begin(), celldofs_a1.end());
      dofs_L.insert(celldofs_L.begin(), celldofs_L.end());

      // Apply zero BC on patch boundary by removing respective dofs
      if (_bc_type == BCType::topological_zero)
      {
        for (FacetIterator facet(*cell); !facet.end(); ++facet)
        {
          if (facet->incident(*vertex))
            continue;

          // FIXME: Factor this code out? facet_dofs are constant across mesh
          std::vector<std::size_t> facet_dofs_a0;
          std::vector<std::size_t> facet_dofs_a1;
          std::vector<std::size_t> facet_dofs_L;
          dofmaps_a[0]->tabulate_facet_dofs(facet_dofs_a0, facet.pos());
          dofmaps_a[1]->tabulate_facet_dofs(facet_dofs_a1, facet.pos());
          (*dofmap_L)[0]->tabulate_facet_dofs(facet_dofs_L, facet.pos());

          for (const auto dof: facet_dofs_a0)
            dofs_a0.erase(celldofs_a0[dof]);
          for (const auto dof: facet_dofs_a1)
            dofs_a1.erase(celldofs_a1[dof]);
          for (const auto dof: facet_dofs_L)
            dofs_L.erase(celldofs_L[dof]);
        }
      }
    }

    // Sort the patch dofs and wrap as ArrayView
    dofs_a0.sort();
    dofs_a1.sort();
    dofs_L.sort();
    dofs_a0_ptr.set(dofs_a0.set());
    dofs_a1_ptr.set(dofs_a1.set());
    dofs_L_ptr.set(dofs_L.set());

    // Check that the local matrix is square
    if (dofs_a0.size() != dofs_a1.size())
    {
      dolfin_error("LocalPatchSolver.cpp",
                   "assemble local LHS",
                   "Local LHS dimensions is non square (%d x %d) on vertex %d",
                   dofs_a0.size(), dofs_a1.size(), vertex->index());
    }

    // Check that the local RHS matches the matrix
    if (dofs_a0.size() != dofs_L.size())
    {
      dolfin_error("LocalPatchSolver.cpp",
                   "assemble local RHS",
                   "Local RHS dimension %d is does not match first dimension "
                   "%d of LHS on vertex %d",
                   dofs_L.size(), dofs_a0.size(), vertex->index());
    }

    // Allocate and zero local tensors (if needed)
    x_e.resize(dofs_a1.size());
    b_e.resize(dofs_L.size(), 1);
    if (!global_b)
      b_e.setZero();
    if (!use_cache)
    {
      A_e.resize(dofs_a0.size(), dofs_a1.size());
      A_e.setZero();
    }

    // Assemble local rhs and/or local matrix cell-by-cell (if needed)
    if (!global_b || !use_cache)
    {
      ufc::cell ufc_cell;
      std::vector<double> coordinate_dofs;

      for (CellIterator cell(*vertex); !cell.end(); ++cell)
      {
        cell->get_coordinate_dofs(coordinate_dofs);

        // Find patch index, i.e., local vertex index in cell
        std::size_t patch;
        {
          VertexIterator v(*cell);
          for (; !v.end(); ++v)
          {
            if (v->index() == vertex->index())
            {
              patch = v.pos();
              break;
            }
          }
          dolfin_assert(!v.end());
        }

        // Assemble local rhs vector if needed
        if (!global_b)
        {
          const ArrayView<const dolfin::la_index> celldofs_L
            = (*dofmap_L)[0]->cell_dofs(cell->index());
          b_cell.resize(celldofs_L.size(), 1);
          LocalAssembler::assemble(b_cell, ufc_L[patch],
                                   coordinate_dofs, ufc_cell, *cell,
                                   cell_domains_L[patch],
                                   exterior_facet_domains_L[patch],
                                   interior_facet_domains_L[patch]);
          const auto _dofmap_L = _compute_cell_to_patch_map(celldofs_L, dofs_L_ptr);
          for (std::size_t i = 0; i < _dofmap_L.size(); i++)
            b_e(_dofmap_L[i], 0) += b_cell(i, 0);
        }

        // Assemble local matrix if needed
        if (!use_cache)
        {
          const ArrayView<const dolfin::la_index> celldofs_a0
            = dofmaps_a[0]->cell_dofs(cell->index());
          const ArrayView<const dolfin::la_index> celldofs_a1
            = dofmaps_a[1]->cell_dofs(cell->index());
          A_cell.resize(celldofs_a0.size(), celldofs_a1.size());
          LocalAssembler::assemble(A_cell, ufc_a[patch],
                                   coordinate_dofs, ufc_cell, *cell,
                                   cell_domains_a[patch],
                                   exterior_facet_domains_a[patch],
                                   interior_facet_domains_a[patch]);
          const auto _dofmap_a0 = _compute_cell_to_patch_map(celldofs_a0, dofs_a0_ptr);
          const auto _dofmap_a1 = _compute_cell_to_patch_map(celldofs_a1, dofs_a1_ptr);
          for (std::size_t i = 0; i < _dofmap_a0.size(); i++)
            for (std::size_t j = 0; j < _dofmap_a1.size(); j++)
              A_e(_dofmap_a0[i], _dofmap_a1[j]) += A_cell(i, j);
        }
      }
    }

    // Copy global RHS data into local RHS vector if needed
    if (global_b)
      // FIXME: This is screwed-up! Why is here index 0?
      (*global_b)[0]->get_local(b_e.data(), dofs_L.size(), dofs_L_ptr.data());

    // Solve local problem
    if (use_cache)
    {
      // Use cached factorisations
      if (_solver_type == SolverType::Cholesky)
        x_e = _cholesky_cache[vertex->index()].solve(b_e);
      else
        x_e = _lu_cache[vertex->index()].solve(b_e);
    }
    else
    {
      // Factorise and solve
      if (_solver_type == SolverType::Cholesky)
      {
        cholesky.compute(A_e);
        x_e = cholesky.solve(b_e);
      }
      else
      {
        lu.compute(A_e);
        x_e = lu.solve(b_e);
      }
    }

    // Add local solution to global vector
    // FIXME: This is screwed-up! Why is here index 0?
    x[0]->add_local(x_e.data(), dofs_a1.size(), dofs_a1_ptr.data());

    // Update progress
    p++;
  }

  // Finalise vector
  // FIXME: This is screwed-up! Why is here index 0?
  x[0]->apply("add");
}
//----------------------------------------------------------------------------
void LocalPatchSolver::factorize()
{
  // FIXME: Implement me
}
//----------------------------------------------------------------------------
void LocalPatchSolver::clear_factorization()
{
  _cholesky_cache.clear();
  _lu_cache.clear();
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::_check_input(
  const std::vector<std::shared_ptr<const Form>>* a,
  const std::vector<std::shared_ptr<const Form>>* L) const
{
  const FunctionSpace* _common_space_0 = nullptr;
  const FunctionSpace* _common_space_1 = nullptr;

  if (a)
  {
    for (const auto ai: *a)
    {
      dolfin_assert(ai);

      // Check rank
      if (ai->rank() != 2)
      {
        dolfin_error("LocalPatchSolver.cpp",
                     "instantiate LocalPatchSolver",
                     "Expected bilinear forms");
      }

      // Check all the forms have same test space
      dolfin_assert(ai->function_space(0));
      if (_common_space_0 && _common_space_0->id() != ai->function_space(0)->id())
      {
        dolfin_error("LocalPatchSolver.cpp",
                     "instantiate LocalPatchSolver",
                     "Expected matching function spaces");
      }
      else
      {
        _common_space_0 = ai->function_space(0).get();
      }

      // Check all the forms have same trial space
      dolfin_assert(ai->function_space(1));
      if (_common_space_1 && _common_space_1->id() != ai->function_space(1)->id())
      {
        dolfin_error("LocalPatchSolver.cpp",
                     "instantiate LocalPatchSolver",
                     "Expected matching function spaces");
      }
      else
      {
        _common_space_1 = ai->function_space(1).get();
      }
    }

    // Check correct number of forms to match patch
    // NOTE: _num_patches has been set before
    dolfin_assert(a->size() == _num_patches);
    if (a->size() == 0 || a->size() != (*a)[0]->mesh()->type().num_vertices())
    {
      dolfin_error("LocalPatchSolver.cpp",
                   "instantiate LocalPatchSolver",
                   "Number of supplied bilinear forms must be equal to number "
                   "number of patches/vertices per cell");
    }
  }

  if (L)
  {
    for (const auto Li: *L)
    {
      dolfin_assert(Li);

      // Check rank
      if (Li->rank() != 1)
      {
        dolfin_error("LocalPatchSolver.cpp",
                     "instantiate LocalPatchSolver",
                     "Expected linear forms");
      }

      // Check all the forms have same space
      dolfin_assert(Li->function_space(0));
      if (_common_space_0 && _common_space_0->id() != Li->function_space(0)->id())
      {
        dolfin_error("LocalPatchSolver.cpp",
                     "instantiate LocalPatchSolver",
                     "Expected matching function spaces");
      }
      else
      {
        _common_space_0 = Li->function_space(0).get();
      }
    }

    // Check correct number of forms to match patch
    // NOTE: _num_patches been set and checked before
    if (L->size() != _num_patches)
    {
      dolfin_error("LocalPatchSolver.cpp",
                   "instantiate LocalPatchSolver",
                   "Number of supplied linear forms must be equal to number "
                   "number of patches/vertices per cell");
    }
  }

  // Check that mesh is common
  dolfin_assert(!_common_space_0 || _common_space_0->mesh());
  dolfin_assert(!_common_space_1 || _common_space_1->mesh());
  if (_common_space_0 && _common_space_1
      && _common_space_0->mesh()->id() != _common_space_1->mesh()->id())
  {
    dolfin_error("LocalPatchSolver.cpp",
                 "instantiate LocalPatchSolver",
                 "Expected common mesh");
  }

  // Check that we are on simplex
  // TODO: Make this simpler by implementing CellType::is_simplex
  const CellType::Type& cell_type = _common_space_0->mesh()->type().cell_type();
  if (!(cell_type == CellType::interval
        || cell_type == CellType::triangle
        || cell_type == CellType::tetrahedron))
  {
    dolfin_error("LocalPatchSolver.cpp",
                 "instantiate LocalPatchSolver",
                 "Expected simplicial mesh");
  }

  // FIXME: Will probably work correctly only on ghosted meshes. Check it?!
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::_check_input(std::vector<Function*> u) const
{
  // Check all the functcion have matching space to trial space
  dolfin_assert(_a.size() > 0 && _a[0] && _a[0]->function_space(1));
  const FunctionSpace& _common_space = *_a[0]->function_space(1);
  for (const auto ui: u)
  {
    dolfin_assert(ui && ui->function_space());
    if (ui->function_space()->id() != _common_space.id())
      dolfin_error("LocalPatchSolver.cpp",
                   "instantiate LocalPatchSolver",
                   "Expected matching function space in function and form");
  }

  // Check correct number of functions to match patch
  // NOTE: _num_patches has been set before
  if (u.size() != _num_patches)
  {
    dolfin_error("LocalPatchSolver.cpp",
                 "use LocalPatchSolver",
                 "Number of supplied functions must be equal to number "
                 "number of patches/vertices per cell");
  }
}
//-----------------------------------------------------------------------------
void LocalPatchSolver::_check_input(std::vector<GenericVector*> x,
                                    std::vector<const GenericVector*> b,
                                    std::vector<const GenericDofMap*> dofmap_b)
  const
{
  // NOTE: _num_patches has been set before
  if (x.size() != _num_patches)
  {
    dolfin_error("LocalPatchSolver.cpp",
                 "use LocalPatchSolver",
                 "Number of supplied solution vectors must be equal to number "
                 "number of patches/vertices per cell");
  }
  if (b.size() != _num_patches)
  {
    dolfin_error("LocalPatchSolver.cpp",
                 "use LocalPatchSolver",
                 "Number of supplied rhs vectors must be equal to number "
                 "number of patches/vertices per cell");
  }
  if (dofmap_b.size() != _num_patches)
  {
    dolfin_error("LocalPatchSolver.cpp",
                 "use LocalPatchSolver",
                 "Number of supplied rhs dofmaps must be equal to number "
                 "number of patches/vertices per cell");
  }

  for (const auto xi: x)
  {
    dolfin_assert(xi);
    // FIXME: Check for compatibility of tensors. We could factor out code
    //        of DirichletBC::check_arguments into DofMap::check_tensor...
  }
  for (const auto bi: b)
  {
    dolfin_assert(bi);
    // FIXME: Check for compatibility of tensors. We could factor out code
    //        of DirichletBC::check_arguments into DofMap::check_tensor...
  }

  // Check all that passed dofmaps match lhs dofmap
  dolfin_assert(_a.size());
  dolfin_assert(_a[0]->function_space(0));
  dolfin_assert(_a[0]->function_space(0)->dofmap());
  const GenericDofMap& _expected_dofmap = *_a[0]->function_space(0)->dofmap();
  for (const auto dofmap_bi: dofmap_b)
  {
    dolfin_assert(dofmap_bi);

    // FIXME: What is the motivation to pass in rhs dofmap in LocalSolver?
    //        Can it be different from lhs dofmap?
    if (dofmap_bi->id() != _expected_dofmap.id())
    {
      dolfin_error("LocalPatchSolver.cpp",
                   "prepare LocalPatchSolver",
                   "Expected matching dofmap on rhs and lhs");
    }
  }
}
//-----------------------------------------------------------------------------
