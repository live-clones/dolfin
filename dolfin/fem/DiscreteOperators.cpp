// Copyright (C) 2015 Garth N. Wells
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

#include <vector>
#include <dolfin/common/ArrayView.h>
#include <dolfin/common/Timer.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/SparsityPatternBuilder.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/TensorLayout.h>
#include <dolfin/log/Progress.h>
#include <dolfin/mesh/Mesh.h>
#include "DiscreteOperators.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::shared_ptr<GenericMatrix>
DiscreteOperators::build_gradient(const FunctionSpace& V0,
                                  const FunctionSpace& V1)
{
  // Get mesh
  dolfin_assert(V0.mesh());
  const Mesh& mesh = *(V0.mesh());

  // Check that mesh is the same for both function spaces
  dolfin_assert(V1.mesh());
  if (&mesh != V1.mesh().get())
  {
    dolfin_error("DiscreteGradient.cpp",
                 "compute discrete gradient operator",
                 "function spaces do not share the same mesh");
  }

  // FIXME: Add checks on finite elements

  // Start timer
  Timer t0("Assemble cells");

  // Declare matrix
  auto A = std::make_shared<Matrix>();

  // Create layout for initialising tensor
  std::shared_ptr<TensorLayout> tensor_layout;
  tensor_layout = A->factory().create_layout(mesh.mpi_comm(), 2);
  dolfin_assert(tensor_layout);

  // Copy index maps from dofmaps
  std::vector<std::shared_ptr<const IndexMap>> index_maps
    = {V0.dofmap()->index_map(), V1.dofmap()->index_map()};
  std::vector<std::pair<std::size_t, std::size_t>> local_range
    = {V0.dofmap()->ownership_range(), V1.dofmap()->ownership_range()};

  // Initialise tensor layout
  tensor_layout->init(index_maps, TensorLayout::Ghosts::UNGHOSTED);

  // Build sparsity pattern (like for cell integral)
  if (tensor_layout->sparsity_pattern())
  {
    SparsityPattern& pattern = *tensor_layout->sparsity_pattern();
    SparsityPatternBuilder::build(pattern, mesh,
                                 {V0.dofmap().get(), V1.dofmap().get()},
                                 true, false, false, false, false);
  }

  // Initialise matrix
  Timer t1("Init matrix");
  A->init(*tensor_layout);
  t1.stop();

  const GenericDofMap& dofmap0 = *V0.dofmap();
  const GenericDofMap& dofmap1 = *V1.dofmap();
  const ufc::finite_element& element0 = *V0.element()->ufc_element();
  const ufc::finite_element& element1 = *V1.element()->ufc_element();

  // Create data structures for local assembly data
  ufc::cell ufc_cell;
  std::vector<double> coordinate_dofs;
  ArrayView<const dolfin::la_index> dofs0;
  ArrayView<const dolfin::la_index> dofs1;
  std::vector<double> values(dofmap0.max_element_dofs()
    * dofmap1.max_element_dofs());
  std::vector<double> column_values(dofmap0.max_element_dofs());
  UFCGradient ufc_gradient(&element1);

  // Assemble over cells
  Progress p("Assembling discrete gradient matrix", mesh.num_cells());
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Check that cell is not a ghost
    dolfin_assert(!cell->is_ghost());

    // Update to current cell
    cell->get_cell_data(ufc_cell);
    cell->get_coordinate_dofs(coordinate_dofs);
    ufc_gradient.update(0, coordinate_dofs.data(), ufc_cell.orientation);

    // Get local-to-global dof maps for cell
    auto dmap0 = dofmap0.cell_dofs(cell->index());
    auto dmap1 = dofmap1.cell_dofs(cell->index());
    dofs0 = ArrayView<const dolfin::la_index>(dmap0.size(), dmap0.data());
    dofs1 = ArrayView<const dolfin::la_index>(dmap1.size(), dmap1.data());

    // Skip if at least one dofmap is empty
    if (dofs0.size() == 0 or dofs1.size() == 0)
      continue;

    // Tabulate dofs on gradients
    for (std::size_t j = 0; j < dofs1.size(); ++j)
    {
      element0.evaluate_dofs(column_values.data(), ufc_gradient,
                             coordinate_dofs.data(), ufc_cell.orientation,
                             ufc_cell);
      for (std::size_t i = 0; i < dofs0.size(); ++i)
        values[j + i*dofs1.size()] = column_values[i];
      ufc_gradient++;
    }

    // Set values in matrix
    A->set_local(values.data(),
                 dofs0.size(), dofs0.data(),
                 dofs1.size(), dofs1.data());

    p++;
  }

  // Finalise matrix
  A->apply("insert");

  return A;
}
//-----------------------------------------------------------------------------
DiscreteOperators::UFCGradient::UFCGradient(const ufc::finite_element* element)
  : _element(element)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
DiscreteOperators::UFCGradient::~UFCGradient()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void DiscreteOperators::UFCGradient::evaluate(double* values,
                                              const double* x,
                                              const ufc::cell& c) const
{
  _element->evaluate_basis_derivatives(_i, 1, values, x, _coordinate_dofs,
                                       _cell_orientation, _cm);
}
//-----------------------------------------------------------------------------
void DiscreteOperators::UFCGradient::update(std::size_t i,
                                            const double* coordinate_dofs,
                                            int cell_orientation,
                                            const ufc::coordinate_mapping* cm)
{
  _i = i;
  _coordinate_dofs = coordinate_dofs;
  _cell_orientation = cell_orientation;
  _cm = cm;
}
//-----------------------------------------------------------------------------
void DiscreteOperators::UFCGradient::operator++(int)
{
  _i++;
}
//-----------------------------------------------------------------------------
