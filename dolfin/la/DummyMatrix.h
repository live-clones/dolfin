// Copyright (C) 2016 Martin Sandve Alnæs
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

#ifndef __DUMMY_MATRIX_H
#define __DUMMY_MATRIX_H

#include <tuple>
#include <vector>
#include "GenericTensor.h"
#include "GenericMatrix.h"
#include "GenericLinearOperator.h"

namespace dolfin
{

  class DummyVector;
  //class TensorLayout;

  /// This class defines a common interface for matrices.

  class DummyMatrix : public GenericMatrix
  {
    MPI_Comm comm;
    std::vector<std::size_t> sizes;

  public:

    /// Create empty matrix (on MPI_COMM_WORLD)
    DummyMatrix();

    /// Create empty matrix
    explicit DummyMatrix(MPI_Comm comm);

    /// Destructor
    virtual ~DummyMatrix();

    //--- Implementation of the GenericLinearOperator interface ---

    /// Compute matrix-vector product y = Ax
    virtual void mult(const GenericVector& x, GenericVector& y) const
    {}

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using tensor layout
    virtual void init(const TensorLayout& tensor_layout);

    /// Return true if empty
    virtual bool empty() const
    { return true; }

    /// Return the MPI communicator
    MPI_Comm mpi_comm() const
    { return comm; }

    /// Return size of given dimension
    virtual std::size_t size(std::size_t dim) const
    { return sizes[dim]; }

    /// Return local ownership range
    virtual std::pair<std::int64_t, std::int64_t>
    local_range(std::size_t dim) const
    { return {0, sizes[dim]}; }

    /// Return number of non-zero entries in matrix (collective)
    virtual std::size_t nnz() const
    { return 0; }

    /// Get block of values
    virtual void get(double* block, const dolfin::la_index* num_rows,
                     const dolfin::la_index * const * rows) const
    {}

    /// Set block of values using global indices
    virtual void set(const double* block, const dolfin::la_index* num_rows,
                     const dolfin::la_index * const * rows)
    {}

    /// Set block of values using local indices
    virtual void set_local(const double* block,
                           const dolfin::la_index* num_rows,
                           const dolfin::la_index * const * rows)
    {}

    /// Add block of values using global indices
    virtual void add(const double* block, const dolfin::la_index* num_rows,
                     const dolfin::la_index * const * rows)
    {}

    /// Add block of values using local indices
    virtual void add_local(const double* block,
                           const dolfin::la_index* num_rows,
                           const dolfin::la_index * const * rows)
    {}

    /// Add block of values using global indices
    virtual void
    add(const double* block,
        const std::vector<ArrayView<const dolfin::la_index>>& rows)
    {}

    /// Add block of values using local indices
    virtual void
    add_local(const double* block,
              const std::vector<ArrayView<const dolfin::la_index>>& rows)
    {}

    /// Set all entries to zero and keep any sparse structure
    virtual void zero()
    {}

    /// Finalize assembly of tensor
    virtual void apply(std::string mode)
    {}

    /// Return informal string representation (pretty-print)
    virtual std::string str(bool verbose) const
    { return "DummyMatrix"; }

    /// Return linear algebra backend factory
    virtual GenericLinearAlgebraFactory& factory() const;

    //--- Matrix interface ---

    /// Return copy of matrix
    virtual std::shared_ptr<GenericMatrix> copy() const
    { return std::shared_ptr<GenericMatrix>(new DummyMatrix(comm)); }

    /// Initialize vector z to be compatible with the matrix-vector
    /// product y = Ax. In the parallel case, both size and layout are
    /// important.
    ///
    /// *Arguments*
    ///     dim (std::size_t)
    ///         The dimension (axis): dim = 0 --> z = y, dim = 1 --> z = x
    virtual void init_vector(GenericVector& z, std::size_t dim) const;

    /// Get block of values
    virtual void get(double* block,
                     std::size_t m, const dolfin::la_index* rows,
                     std::size_t n, const dolfin::la_index* cols) const
    {}

    /// Set block of values using global indices
    virtual void set(const double* block,
                     std::size_t m, const dolfin::la_index* rows,
                     std::size_t n, const dolfin::la_index* cols)
    {}

    /// Set block of values using local indices
    virtual void set_local(const double* block,
                           std::size_t m, const dolfin::la_index* rows,
                           std::size_t n, const dolfin::la_index* cols)
    {}

    /// Add block of values using global indices
    virtual void add(const double* block,
                     std::size_t m, const dolfin::la_index* rows,
                     std::size_t n, const dolfin::la_index* cols)
    {}

    /// Add block of values using local indices
    virtual void add_local(const double* block,
                           std::size_t m, const dolfin::la_index* rows,
                           std::size_t n, const dolfin::la_index* cols)
    {}

    /// Add multiple of given matrix (AXPY operation)
    virtual void axpy(double a, const GenericMatrix& A,
                      bool same_nonzero_pattern)
    {}

    /// Return norm of matrix
    virtual double norm(std::string norm_type) const
    { return 0.0; }

    /// Get non-zero values of given row (global index) on local process
    virtual void getrow(std::size_t row, std::vector<std::size_t>& columns,
                        std::vector<double>& values) const
    {}

    /// Set values for given row (global index) on local process
    virtual void setrow(std::size_t row,
                        const std::vector<std::size_t>& columns,
                        const std::vector<double>& values)
    {}

    /// Set given rows (global row indices) to zero
    virtual void zero(std::size_t m, const dolfin::la_index* rows)
    {}

    /// Set given rows (local row indices) to zero
    virtual void zero_local(std::size_t m, const dolfin::la_index* rows)
    {}

    /// Set given rows (global row indices) to identity matrix
    virtual void ident(std::size_t m, const dolfin::la_index* rows)
    {}

    /// Set given rows (local row indices) to identity matrix
    virtual void ident_local(std::size_t m, const dolfin::la_index* rows)
    {}

    /// Matrix-vector product, y = A^T x. The y vector must either be
    /// zero-sized or have correct size and parallel layout.
    virtual void transpmult(const GenericVector& x, GenericVector& y) const
    {}

    /// Get diagonal of a matrix
    virtual void get_diagonal(GenericVector& x) const
    {}

    /// Set diagonal of a matrix
    virtual void set_diagonal(const GenericVector& x)
    {}

    /// Multiply matrix by given number
    virtual const GenericMatrix& operator*= (double a)
    { return *this; }

    /// Divide matrix by given number
    virtual const GenericMatrix& operator/= (double a)
    { return *this; }

    /// Add given matrix
    const GenericMatrix& operator+= (const GenericMatrix& A)
    { return *this; }

    /// Subtract given matrix
    const GenericMatrix& operator-= (const GenericMatrix& A)
    { return *this; }

    /// Test if matrix is symmetric
    virtual bool is_symmetric(double tol) const
    {
      dolfin_error("DummyMatrix.h",
                   "test if matrix is symmetric",
                   "Not implemented by current linear algebra backend");
      return false;
    }

    /// Assignment operator
    virtual const GenericMatrix& operator= (const GenericMatrix& x)
    { return *this; }

    //--- Convenience functions ---

    /// Get value of given entry
    virtual double operator() (dolfin::la_index i, dolfin::la_index j) const
    { return 0.0; }

    /// Get value of given entry
    virtual double getitem(std::pair<dolfin::la_index,
                           dolfin::la_index> ij) const
    { return 0.0; }

    /// Set given entry to value. apply("insert") must be called
    /// before using using the object.
    virtual void setitem(std::pair<dolfin::la_index, dolfin::la_index> ij,
                         double value)
    {}

    /// Insert one on the diagonal for all zero rows
    virtual void ident_zeros()
    {}

  };

}

#endif
