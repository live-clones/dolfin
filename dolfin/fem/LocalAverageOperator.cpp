// Copyright (C) 2013 Mikael Mortensen
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
// First added:
// Last changed:

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Timer.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/Vector.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include "GenericDofMap.h"
#include "Form.h"
#include "UFC.h"

#include "LocalAverageOperator.h"

using namespace dolfin;

//----------------------------------------------------------------------------
LocalAverageOperator::LocalAverageOperator(const FunctionSpace& V) :
  weight(new Function(V))
{
  compute_avg_weight(V);
}

LocalAverageOperator::LocalAverageOperator(const FunctionSpace& V, const Form& test) :
  weight(new Function(V)), test(reference_to_no_delete_pointer(test))
{
  compute_lumped_weight(V);
}

boost::shared_ptr<Function> LocalAverageOperator::get_weight()
{
  return weight;
}

void LocalAverageOperator::compute_avg_weight(const FunctionSpace& V)
{
  // Compute weights for averaging with neighboring cells
    
  // Get the mesh, element and dofmap
  boost::shared_ptr<const Mesh> mesh = V.mesh();
  boost::shared_ptr<const FiniteElement> element = V.element();
  boost::shared_ptr<const GenericDofMap> dofmap_u = V.dofmap();

  // Allocate storage for weights on one cell
  std::vector<double> ws(element->space_dimension()); 
  
  // Compute weights
  GenericVector& w_vector = *weight->vector();  
  for (CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    const std::vector<dolfin::la_index>& dofs
      = dofmap_u->cell_dofs(cell->index());
      
    std::fill(ws.begin(), ws.end(), 1./cell->volume());
//    std::fill(ws.begin(), ws.end(), 1.);
    w_vector.add(ws.data(), dofs.size(), dofs.data());      
  }  
  w_vector.apply("insert");  
  
  // Get the ownership range
  std::pair<std::size_t, std::size_t> local_range = 
    dofmap_u->ownership_range();  

  // Invert the weights 
  std::vector<double> local_x;
  w_vector.get_local(local_x);    
  for (std::size_t i = 0; i < local_x.size(); i++)            
      w_vector.setitem(i+local_range.first, 1./local_x[i]);
}

void LocalAverageOperator::compute_lumped_weight(const FunctionSpace& V)
{
  // Compute weights for a lumped average with neighboring cells
  UFC ufc_b(*test);
 
  // Get the mesh, element and dofmap
  boost::shared_ptr<const Mesh> mesh = V.mesh();
  boost::shared_ptr<const FiniteElement> element = V.element();
  boost::shared_ptr<const GenericDofMap> dofmap_u = V.dofmap();

  ufc::cell_integral* integral_b = ufc_b.default_cell_integral.get();
  
  // Compute weights
  GenericVector& w_vector = *weight->vector();  
  for (CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    ufc_b.update(*cell);

    const std::vector<dolfin::la_index>& dofs
      = dofmap_u->cell_dofs(cell->index());
      
    integral_b->tabulate_tensor(&ufc_b.A[0],
                                ufc_b.w(),
                                &ufc_b.cell.vertex_coordinates[0],
                                ufc_b.cell.orientation);
    
    w_vector.add(&ufc_b.A[0], dofs.size(), dofs.data());
    
  }
  w_vector.apply("insert");  
  
  // Get the ownership range
  std::pair<std::size_t, std::size_t> local_range = 
    dofmap_u->ownership_range();  

  // Invert the weights 
  std::vector<double> local_x;
  w_vector.get_local(local_x);    
  for (std::size_t i = 0; i < local_x.size(); i++)            
      w_vector.setitem(i+local_range.first, 1./local_x[i]);  
}

void LocalAverageOperator::solve(Function& u, const Form& a) const
{
  // This version works for scalar function spaces only
  // Example usage:
  //    V = FunctionSpace(mesh, 'CG', 1)
  //    u = Function(V)
  //    du= Function(V)
  //    lp = LocalAverageOperator(V)
  //    lp.solve(du, u.dx(0)*dx())
  //  
  assert(u.function_space() == weight->function_space());
  
  UFC ufc_a(a);  

  // Set timer
  Timer timer("Local average operator");

  // Extract mesh
  const Mesh& mesh = a.mesh();
  
  // Update off-process coefficients
  const std::vector<boost::shared_ptr<const GenericFunction> >
    coefficients_a = a.coefficients();
  for (std::size_t i = 0; i < coefficients_a.size(); ++i)
    coefficients_a[i]->update();

  // Form ranks
  const std::size_t rank_a = ufc_a.form.rank();

  // Check form ranks
  dolfin_assert(rank_a == 0);
  
    // Get element and dofmap
  boost::shared_ptr<const FiniteElement> element 
    = u.function_space()->element();
  boost::shared_ptr<const GenericDofMap> dofmap_u 
    = u.function_space()->dofmap();

  // Cell integrals
  ufc::cell_integral* integral_a = ufc_a.default_cell_integral.get();
    
  // Allocate vector to store dofs on cell
  std::vector<double> vals(element->space_dimension());
  
  // Get vector pointer
  GenericVector& u_vector = *u.vector();
  
  // Assemble over cells
  Progress p("Performing local averaging", mesh.num_cells());
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_a.update(*cell);
    
    const std::vector<dolfin::la_index>& dofs
      = dofmap_u->cell_dofs(cell->index());

    // Tabulate functional on cell
    integral_a->tabulate_tensor(&ufc_a.A[0],
                                ufc_a.w(),
                                &ufc_a.cell.vertex_coordinates[0],
                                ufc_a.cell.orientation);
        
    std::fill(vals.begin(), vals.end(), ufc_a.A[0]/(pow(cell->volume(), 2)));
    u_vector.add(vals.data(), dofs.size(), dofs.data());
    
    p++;
  }
  u_vector.apply("insert");  
  
  // Scale tensor by weights  
  u_vector *= (*weight->vector());
}

//------------------------------------------------------------------//
void LocalAverageOperator::solve_vector(Function& u, const Form& a) const
{
  // This version is supposed to work also for vector function spaces.
  // Example usage:
  //    V = FunctionSpace(mesh, 'CG', 1)
  //    VV= VectorFunctionSpace(mesh, 'CG', 1)
  //    R = VectorFunctionSpace(mesh, 'R', 0)
  //    c = TestFunction(R)
  //    u = Function(V)
  //    du= Function(VV)
  //    lp = LocalAverageOperator(V)
  //    lp.solve(du, dot(grad(u), c)*dx())
  //
  // Should also work for
  //    du= Function(V)
  //    lp = LocalAverageOperator(V)
  //    lp.solve(du, u.dx(0)*dx())
  //
  UFC ufc_a(a);

  // Set timer
  Timer timer("Local average operator vector");

  // Extract mesh
  const Mesh& mesh = a.mesh();
  
  // Update off-process coefficients
  const std::vector<boost::shared_ptr<const GenericFunction> >
    coefficients_a = a.coefficients();
  for (std::size_t i = 0; i < coefficients_a.size(); ++i)
    coefficients_a[i]->update();
  
  // Get function space, element and dofmap
  boost::shared_ptr<const FunctionSpace> V = u.function_space();
  boost::shared_ptr<const FiniteElement> element = V->element();
  boost::shared_ptr<const GenericDofMap> dofmap_u = V->dofmap();

  // Cell integrals
  ufc::cell_integral* integral_a = ufc_a.default_cell_integral.get();
  
  std::size_t value_size_loc = 1;
  for (std::size_t i = 0; i < element->value_rank(); i++)
     value_size_loc *= element->value_dimension(i);

  std::vector<double> vals(element->space_dimension());
  std::vector<double> ws(element->space_dimension());

  GenericVector& u_vector = *u.vector();  
  // Assemble over cells
  Progress p("Performing local averaging", mesh.num_cells());
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_a.update(*cell);
    
    const std::vector<dolfin::la_index>& dofs
      = dofmap_u->cell_dofs(cell->index());

    // Tabulate functional on cell
    integral_a->tabulate_tensor(&ufc_a.A[0],
                                ufc_a.w(),
                                &ufc_a.cell.vertex_coordinates[0],
                                ufc_a.cell.orientation);
        
    double cv = cell->volume();
    std::vector<double>::iterator it = vals.begin();
    if (element->num_sub_elements() > 0)
    {
      for (std::size_t i = 0; i < value_size_loc; i++) 
      {
        std::size_t end = (*V)[i]->element()->space_dimension();
        std::fill(it, it+end, ufc_a.A[i]/(cv*cv));
        it += end;
      }
    }
    else
    {
      std::fill(vals.begin(), vals.end(), ufc_a.A[0]/(cv*cv));
    }
    
    u_vector.add(vals.data(), dofs.size(), dofs.data());
    
    p++;
  }  
  u_vector.apply("insert");  
  
  // Scale tensor by weight  
  u_vector *= (*weight->vector());
}
//-----------------------------------------------------------------------------

void LocalAverageOperator::solve_lumping(Function& u, const Form& a) const
{
  UFC ufc_a(a);
  
  // Set timer
  Timer timer("Local average operator lumping");

  // Extract mesh
  const Mesh& mesh = a.mesh();
  
  // Update off-process coefficients
  const std::vector<boost::shared_ptr<const GenericFunction> >
    coefficients_a = a.coefficients();
  for (std::size_t i = 0; i < coefficients_a.size(); ++i)
    coefficients_a[i]->update();

  // Form ranks
  const std::size_t rank_a = ufc_a.form.rank();

  // Check form ranks
  dolfin_assert(rank_a == 1);
  
  // Get function space, element and dofmap
  boost::shared_ptr<const FunctionSpace> V = u.function_space();
  boost::shared_ptr<const FiniteElement> element = V->element();
  boost::shared_ptr<const GenericDofMap> dofmap_u = V->dofmap();

  // Cell integrals
  ufc::cell_integral* integral_a = ufc_a.default_cell_integral.get();
  
  std::vector<double> vals(element->space_dimension());
  std::vector<double> ws(element->space_dimension());
  
  GenericVector& u_vector = *u.vector();
  
  // Assemble over cells
  Progress p("Performing local averaging", mesh.num_cells());
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_a.update(*cell);
    
    const std::vector<dolfin::la_index>& dofs
      = dofmap_u->cell_dofs(cell->index());

    // Tabulate functional on cell
    integral_a->tabulate_tensor(&ufc_a.A[0],
                                ufc_a.w(),
                                &ufc_a.cell.vertex_coordinates[0],
                                ufc_a.cell.orientation);
       
    u_vector.add(&ufc_a.A[0], dofs.size(), dofs.data());
    
    p++;
  }  
  u_vector.apply("insert");  

  // Scale tensor by weight  
  u_vector *= (*weight->vector());  
}
