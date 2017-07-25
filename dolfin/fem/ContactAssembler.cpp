
#include <memory>

#include <dolfin/common/Timer.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/GenericLinearAlgebraFactory.h>
#include <dolfin/la/TensorLayout.h>
#include <dolfin/log/log.h>
#include <dolfin/common/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>

#include "FiniteElement.h"
#include "Form.h"
#include "GenericDofMap.h"
#include "SparsityPatternBuilder.h"
#include "AssemblerBase.h"

#include "ContactAssembler.h"



using namespace dolfin;

//-----------------------------------------------------------------------------
ContactAssembler::ContactAssembler(
    const GeometricContact& gc)
: _gc(gc)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ContactAssembler::assemble(GenericTensor& A, const Form& a)
{

  // Initialize global tensor
  _init_global_tensor(A, a);

  Assembler::assemble(A, a);
}

//-----------------------------------------------------------------------------
void ContactAssembler::_init_global_tensor(GenericTensor& A, const Form& a)
{
  dolfin_assert(a.ufc_form());

  // Get dof maps
  std::vector<const GenericDofMap*> dofmaps;
  for (std::size_t i = 0; i < a.rank(); ++i)
    dofmaps.push_back(a.function_space(i)->dofmap().get());

  if (A.empty())
  {
    Timer t0("Build sparsity");

    // Create layout for initialising tensor
    std::shared_ptr<TensorLayout> tensor_layout;
    tensor_layout = A.factory().create_layout(a.rank());
    dolfin_assert(tensor_layout);

    // Get dimensions and mapping across processes for each dimension
    std::vector<std::shared_ptr<const IndexMap>> index_maps;
    for (std::size_t i = 0; i < a.rank(); i++)
    {
      dolfin_assert(dofmaps[i]);
      index_maps.push_back(dofmaps[i]->index_map());
    }

    // Get mesh
    dolfin_assert(a.mesh());
    const Mesh& mesh = *(a.mesh());

    // Initialise tensor layout
    // FIXME: somewhere need to check block sizes are same on both axes
    // NOTE: Jan: that will be done on the backend side; IndexMap will
    //            provide tabulate functions with arbitrary block size;
    //            moreover the functions will tabulate directly using a
    //            correct int type

    tensor_layout->init(mesh.mpi_comm(), index_maps,
                        TensorLayout::Ghosts::UNGHOSTED);

    // Build sparsity pattern if required
    if (tensor_layout->sparsity_pattern())
    {
      SparsityPattern& pattern = *tensor_layout->sparsity_pattern();

      // Build sparsity pattern for problem without contact initially
      SparsityPatternBuilder::build(pattern, mesh,
                                    dofmaps,
                                    a.ufc_form()->has_cell_integrals(),
                                    a.ufc_form()->has_interior_facet_integrals(),
                                    a.ufc_form()->has_exterior_facet_integrals(),
                                    a.ufc_form()->has_vertex_integrals(),
                                    keep_diagonal,
                                    true, false);

      // Additional contact pair sparsity
      SparsityPatternBuilder::build_contact_sparsity_pattern(pattern,
                                    mesh, dofmaps, _gc,
                                    a.ufc_form()->has_cell_integrals(),
                                    a.ufc_form()->has_interior_facet_integrals(),
                                    a.ufc_form()->has_exterior_facet_integrals(),
                                    a.ufc_form()->has_vertex_integrals(),
                                    keep_diagonal);
    }
    t0.stop();

    // Initialize tensor
    Timer t1("Init tensor");
    A.init(*tensor_layout);
    t1.stop();

    // Insert zeros on the diagonal as diagonal entries may be
    // prematurely optimised away by the linear algebra backend when
    // calling GenericMatrix::apply, e.g. PETSc does this then errors
    // when matrices have no diagonal entry inserted.
    if (A.rank() == 2 && keep_diagonal)
    {
      // Down cast to GenericMatrix
      GenericMatrix& _matA = A.down_cast<GenericMatrix>();

      // Loop over rows and insert 0.0 on the diagonal
      const double block = 0.0;
      const std::pair<std::size_t, std::size_t> row_range = A.local_range(0);
      const std::size_t range = std::min(row_range.second, A.size(1));
      for (std::size_t i = row_range.first; i < range; i++)
      {
        dolfin::la_index _i = i;
        _matA.set(&block, 1, &_i, 1, &_i);
      }
      A.apply("flush");
    }

    // Delete sparsity pattern
    Timer t2("Delete sparsity");
    t2.stop();
  }
  else
  {
    // If tensor is not reset, check that dimensions are correct
    for (std::size_t i = 0; i < a.rank(); ++i)
    {
      if (A.size(i) != dofmaps[i]->global_dimension())
      {
        dolfin_error("AssemblerBase.cpp",
                     "assemble form",
                     "Dim %d of tensor does not match form", i);
      }
    }
  }

  if (!add_values)
    A.zero();
}