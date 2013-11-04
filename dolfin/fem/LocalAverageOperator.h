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

#ifndef __LOCAL_AVERAGE_OPERATOR_H
#define __LOCAL_AVERAGE_OPERATOR_H

namespace dolfin
{
    
  /// This class contains routines for performing "cheap" projections.
  /// The idea is to assemble and solve locally by averaging over 
  /// neighbouring cells, without going through a global coefficient     
  /// matrix and a global solve. 
  ///
  /// The accuracy is probably best for projecting piecewise constant or
  /// linear fields into linear Lagrange elements. For higher order 
  /// elements it is probably not accurate enough, but it will still work.
  ///
  /// The algorithm takes a zero rank Form as argument. 
  ///    1) The zero rank Form is assembled for all cells and the constant
  ///       value for cell i is ci and the volume of the cell is wi.
  ///       The value ci/wi is the approximation of the form on the cell.
  ///    2) Add ci/wi * (1/wi) to each dof on the cell for the FunctionSpace
  ///       we are "projecting" into. (1/wi) is the weight chosen such
  ///       that small cells receive large weights.
  ///       Add also (1/wi) to a global weight (Coefficient) on the same
  ///       FunctionSpace as the solution.
  ///    3) After running through all cells, divide the solution by the 
  ///       total weight.
  /// 
  /// For uniform meshes the algorithm is identical to computing dof 
  /// values using the arithmetic average over all cells sharing the dof.
  ///
  /// Usage to get p.dx(0) onto linear Lagrange elements:
  ///   V = FunctionSpace(mesh, 'CG', 1)
  ///   lp = LocalAverageOperator(V)
  ///   p = interpolate(...) # Some solution
  ///   dpdx = Function(V)
  ///   lp.solve(p, p.dx(0)*dx())
  ///
  /// The solution will be very similar to dpdx = project(p.dx(0), V)
  ///
    
  // Forward declarations
  class Function;
  class Form;

  class LocalAverageOperator
  {
  public:
      
    /// Constructors
    LocalAverageOperator(const FunctionSpace& V);

    LocalAverageOperator(const FunctionSpace& V, const Form& test);

    /// Solve by performing local averaging of Form a 
    void solve(Function& u, const Form& a) const;

    /// Solve possible for vector function spaces as well.
    void solve_vector(Function& u, const Form& a) const;
    
    /// Solve using lumped weights. That is, we use the test function
    /// as weight and assemble(TestFunction(V)*dx) as weights.
    /// Solution is computed using 
    ///   v = TestFunction(V)
    ///   lp = LocalAverageOperator(V, v*dx)
    ///   lp.solve(p, v*p.dx(0)*dx())
    void solve_lumping(Function& u, const Form& a) const;
    
    void solve_ls(GenericVector& x, const Form& a, const Form& L, bool symmetric=false) const;
        
    /// Weights are precomputed on construction
    void compute_avg_weight(const FunctionSpace& V);

    void compute_lumped_weight(const FunctionSpace& V);
    
    boost::shared_ptr<Function> get_weight();

  private:
      
    /// Storage for the precomputed weights
    boost::shared_ptr<Function> weight;
    
    /// test Form should be v*dx
    boost::shared_ptr<const Form> test;
    
  };

}

#endif
