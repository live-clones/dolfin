// Copyright (C) 2008-2013 Anders Logg and Garth N. Wells
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

#ifndef __FINITE_ELEMENT_H
#define __FINITE_ELEMENT_H

#include <cmath>
#include <memory>
#include <vector>
#include <ufc.h>
#include <boost/multi_array.hpp>
#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace dolfin
{

  class Cell;

  /// This is a wrapper for a UFC finite element (ufc::finite_element).

  class FiniteElement
  {
  public:

    /// Create finite element from UFC finite element (data may be shared)
    /// @param element (ufc::finite_element)
    ///  UFC finite element
    FiniteElement(std::shared_ptr<const ufc::finite_element> element);

    /// Destructor
    virtual ~FiniteElement() {}

    //--- Direct wrappers for ufc::finite_element ---

    /// Return a string identifying the finite element
    /// @return std::string
    std::string signature() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->signature();
    }

    /// Return the cell shape
    /// @return ufc::shape
    ufc::shape cell_shape() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->cell_shape();
    }

    /// Return the topological dimension of the cell shape
    /// @return std::size_t
    std::size_t topological_dimension() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->topological_dimension();
    }

    /// Return the geometric dimension of the cell shape
    /// @return unsigned int
    virtual unsigned int geometric_dimension() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->geometric_dimension();
    }

    /// Return the dimension of the finite element function space
    /// @return std::size_t
    std::size_t space_dimension() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->space_dimension();
    }

    /// Return the rank of the value space
    std::size_t value_rank() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->value_rank();
    }

    /// Return the dimension of the value space for axis i
    std::size_t value_dimension(std::size_t i) const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->value_dimension(i);
    }

    /// Evaluate basis functions on reference element
    ///
    /// @param[in,out] values (boost::multi_array<double, 3>)
    ///         The basis function values at the points. Shape is
    ///         [num_points][num_dofs][reference_value_size]. Will be
    ///         resized, if required.
    /// @param[in]    X (boost::multi_array<double, 2>&)
    ///         Points on the reference cell at which basis functions are
    ///         evaluated. Shape is [num_points][topological dim]
    void evaluate_reference_basis(boost::multi_array<double, 3>& values,
                                  const boost::multi_array<double, 2>& X) const
    {
      dolfin_assert(_ufc_element);

      // Check element topological dimension is consistent with point
      // type
      assert(X.shape()[1] == _ufc_element->topological_dimension());

      // Get number of points
      const std::size_t num_points = X.shape()[0];

      // Resize values
      std::size_t space_dim = _ufc_element->space_dimension();
      std::size_t reference_value_size = _ufc_element->reference_value_size();
      values.resize(boost::extents[num_points][space_dim][reference_value_size]);

      // Tabulate
      _ufc_element->evaluate_reference_basis(values.data(), num_points, X.data());
    }

    /// Evaluate basis functions on reference element
    ///
    /// @param[in,out] values (boost::multi_array<double, 4>&)
    ///         The basis function values at the points. Shape is
    ///         [num_points][num_dofs][num_derivs][reference_value_size]. Will
    ///         be resized, if required.
    /// @param[in]    order (std::size_t)
    ///         Derivative order
    /// @param[in]    X (boost::multi_array<double, 2>&)
    ///         Points on the reference cell at which basis functions are
    ///         evaluated. Shape is [num_points][topological dim]
    void evaluate_reference_basis_derivatives(boost::multi_array<double, 4>& values,
                                              std::size_t order,
                                              const boost::multi_array<double, 2>& X) const
    {
      dolfin_assert(_ufc_element);

      // Check element topological dimension is consistent with point
      // type
      std::size_t tdim = _ufc_element->topological_dimension();
      assert(X.shape()[1] == tdim);

      // Get number of points
      const std::size_t num_points = X.shape()[0];

      // Number of derivatives
      const std::size_t num_derivs = std::pow(tdim, order);

      // Resize values
      std::size_t space_dim = _ufc_element->space_dimension();
      std::size_t reference_value_size = _ufc_element->reference_value_size();
      values.resize(boost::extents[num_points][space_dim][num_derivs][reference_value_size]);

      // Tabulate
      _ufc_element->evaluate_reference_basis_derivatives(values.data(), order,
                                                         num_points, X.data());
    }


    /// Transform basis functions (derivatives) on reference element
    /// to physical space (push forward). Use degree 0 for Piola
    /// transformations for bases.
    ///
    /// @param[in,out] values (boost::multi_array<double, 4>)
    ///         The transformed values. Shape is
    ///         [num_points][num_dofs][num_derivs][value_size]. Will
    ///         be resized, if required.
    /// @param[in]    order (std::size_t)
    ///         Derivative order
    /// @param[in] reference_values (boost::multi_array<double, 4>)
    ///         Values on the reference element. Shape is
    ///         [num_points][num_dofs][num_derivs][reference_value_size]. Will
    /// @param[in]    X (boost::multi_array<double, 2>&)
    ///         Points on the reference cell at which basis functions are
    ///         evaluated. Shape is [num_points][topological dim]
    /// @param[in]    J (boost::multi_array<double, 3>&)
    ///         Jacobian of the transformation, dx/dX.
    ///         Shape is [num_points][geometric dim][topological dim]
    /// @param[in]    detJ (boost::multi_array<double, 1>&)
    ///         Determinant of the Jacobian. Shape is [num_points].
    /// @param[in]    K (boost::multi_array<double, 3>&)
    ///         (Pseudo)-Inverse of the Jacobian
    ///         Shape is [num_points][topological dim][geometric dim]
    /// @param[in]    cell_orientation (unsigned int)
    ///         Orientation of the cell, 1 means flipped w.r.t. reference cell.
    ///         Only relevant on manifolds (tdim < gdim).
    void transform_reference_basis_derivatives(boost::multi_array<double, 4>& values,
                                               std::size_t order,
                                               const boost::multi_array<double, 4>& reference_values,
                                               const boost::multi_array<double, 2>& X,
                                               const boost::multi_array<double, 3>& J,
                                               const boost::multi_array<double, 1>& detJ,
                                               const boost::multi_array<double, 3>& K,
                                               int cell_orientation)
    {
      dolfin_assert(_ufc_element);

      const std::size_t tdim = _ufc_element->topological_dimension();
      const std::size_t gdim = _ufc_element->geometric_dimension();
      const std::size_t space_dim = _ufc_element->space_dimension();
      const std::size_t value_size = _ufc_element->value_size();
      const std::size_t ref_value_size = _ufc_element->reference_value_size();

      // Number of derivatives
      const std::size_t num_derivs = std::pow(tdim, order);

      // Number of points
      const std::size_t num_points = X.shape()[0];

      // Check dimesions
      dolfin_assert(values.shape()[0] == num_points);
      dolfin_assert(values.shape()[1] == space_dim);
      dolfin_assert(values.shape()[2] == num_derivs);  // FIXME: check if to should be gdim**order or tdim**order
      dolfin_assert(values.shape()[3] == value_size);

      dolfin_assert(reference_values.shape()[0] == num_points);
      dolfin_assert(reference_values.shape()[1] == space_dim);
      dolfin_assert(reference_values.shape()[2] == num_derivs);
      dolfin_assert(reference_values.shape()[3] == ref_value_size);

      dolfin_assert(X.shape()[1] == tdim);

      dolfin_assert(J.shape()[0] == num_points);
      dolfin_assert(J.shape()[1] == gdim);
      dolfin_assert(J.shape()[2] == tdim);

      dolfin_assert(detJ.shape()[0] == num_points);

      dolfin_assert(K.shape()[0] == num_points);
      dolfin_assert(K.shape()[1] == tdim);
      dolfin_assert(K.shape()[2] == gdim);

      // Resize values array
      values.resize(boost::extents[num_points][space_dim][num_derivs][value_size]);

      // Transform data
      _ufc_element->transform_reference_basis_derivatives(values.data(),
                                                          order,
                                                          num_points,
                                                          reference_values.data(),
                                                          X.data(),
                                                          J.data(),
                                                          detJ.data(),
                                                          K.data(),
                                                          cell_orientation);
    }

    /// Evaluate basis function i at given point in cell
    void evaluate_basis(std::size_t i, double* values, const double* x,
                        const double* coordinate_dofs,
                        int cell_orientation) const
    {
      dolfin_assert(_ufc_element);
      _ufc_element->evaluate_basis(i, values, x, coordinate_dofs,
                                   cell_orientation);
    }

    /// Evaluate all basis functions at given point in cell
    void evaluate_basis_all(double* values,
                            const double* x,
                            const double* coordinate_dofs,
                            int cell_orientation) const
    {
      dolfin_assert(_ufc_element);
      _ufc_element->evaluate_basis_all(values, x, coordinate_dofs,
                                       cell_orientation);
    }

    // DEPRECATED
    /// Evaluate order n derivatives of basis function i at given
    /// point in cell
    void evaluate_basis_derivatives(unsigned int i,
                                    unsigned int n,
                                    double* values,
                                    const double* x,
                                    const double* coordinate_dofs,
                                    int cell_orientation) const
    {
      dolfin_assert(_ufc_element);
      _ufc_element->evaluate_basis_derivatives(i, n, values, x,
                                               coordinate_dofs,
                                               cell_orientation);
    }

    // DEPRECATED
    /// Evaluate order n derivatives of all basis functions at given
    /// point in cell
    void evaluate_basis_derivatives_all(unsigned int n,
                                        double* values,
                                        const double* x,
                                        const double* coordinate_dofs,
                                        int cell_orientation) const
    {
      dolfin_assert(_ufc_element);
      _ufc_element->evaluate_basis_derivatives_all(n, values, x,
                                                   coordinate_dofs,
                                                   cell_orientation);
    }

    /// Evaluate linear functional for dof i on the function f
    double evaluate_dof(std::size_t i,
                        const ufc::function& function,
                        const double* coordinate_dofs,
                        int cell_orientation,
                        const ufc::cell& c) const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->evaluate_dof(i, function, coordinate_dofs,
                                        cell_orientation, c);
    }

    /// Evaluate linear functionals for all dofs on the function f
    void evaluate_dofs(double* values,
                       const ufc::function& f,
                       const double* coordinate_dofs,
                       int cell_orientation,
                       const ufc::cell& c) const
    {
      dolfin_assert(_ufc_element);
      _ufc_element->evaluate_dofs(values, f, coordinate_dofs,
                                  cell_orientation, c);
    }

    /// Interpolate vertex values from dof values
    /// @param vertex_values (double*)
    /// @param coefficients (double*)
    /// @param coordinate_dofs (const double*)
    /// @param cell_orientation (int)
    void interpolate_vertex_values(double* vertex_values,
                                   double* coefficients,
                                   const double* coordinate_dofs,
                                   int cell_orientation) const
    {
      dolfin_assert(_ufc_element);
      _ufc_element->interpolate_vertex_values(vertex_values, coefficients,
                                              coordinate_dofs,
                                              cell_orientation);
    }

    /// Tabulate the coordinates of all dofs on an element
    ///
    /// @param[in,out]    coordinates (boost::multi_array<double, 2>)
    ///         The coordinates of all dofs on a cell.
    /// @param[in]    coordinate_dofs (std::vector<double>)
    ///         The cell coordinates
    /// @param[in]    cell (Cell)
    ///         The cell.
    void tabulate_dof_coordinates(boost::multi_array<double, 2>& coordinates,
                                  const std::vector<double>& coordinate_dofs,
                                  const Cell& cell) const;

    /// Return the number of sub elements (for a mixed element)
    /// @return std::size_t
    ///   number of sub-elements
    std::size_t num_sub_elements() const
    {
      dolfin_assert(_ufc_element);
      return _ufc_element->num_sub_elements();
    }

    //--- DOLFIN-specific extensions of the interface ---

    /// Return simple hash of the signature string
    std::size_t hash() const
    { return _hash; }

    /// Create a new finite element for sub element i (for a mixed element)
    std::shared_ptr<const FiniteElement>
      create_sub_element(std::size_t i) const
    {
      dolfin_assert(_ufc_element);
      std::shared_ptr<const ufc::finite_element>
        ufc_element(_ufc_element->create_sub_element(i));
      std::shared_ptr<const FiniteElement>
        element(new const FiniteElement(ufc_element));
      return element;
    }

    /// Create a new class instance
    std::shared_ptr<const FiniteElement> create() const
    {
      dolfin_assert(_ufc_element);
      std::shared_ptr<const ufc::finite_element>
        ufc_element(_ufc_element->create());
      return std::shared_ptr<const FiniteElement>(new FiniteElement(ufc_element));
    }

    /// Extract sub finite element for component
    std::shared_ptr<const FiniteElement>
      extract_sub_element(const std::vector<std::size_t>& component) const;

    /// Return underlying UFC element. Intended for libray usage only
    /// and may change.
    std::shared_ptr<const ufc::finite_element> ufc_element() const
    { return _ufc_element; }

  private:

    // UFC finite element
    std::shared_ptr<const ufc::finite_element> _ufc_element;

    // Recursively extract sub finite element
    static std::shared_ptr<const FiniteElement>
      extract_sub_element(const FiniteElement& finite_element,
                          const std::vector<std::size_t>& component);

    // Simple hash of the signature string
    std::size_t _hash;

  };

}
#endif
