// Copyright (C) 2016 Chris Richardson
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

#ifndef __PETSC_NESTMATRIX_H
#define __PETSC_NESTMATRIX_H

#ifdef HAS_PETSC

#include <petscmat.h>
#include <petscsys.h>

#include "PETScMatrix.h"

namespace dolfin
{
  class GenericMatrix;
  class GenericVector;

  class PETScNestMatrix : public PETScMatrix
  {
  public:

    /// Create empty matrix
    PETScNestMatrix();

    /// Create from a list of matrices
    explicit PETScNestMatrix
      (std::vector<std::shared_ptr<const GenericMatrix>> mats);

    /// Destructor
    virtual ~PETScNestMatrix();

    /// Multiply
    virtual void mult(const GenericVector& x, GenericVector& y) const;

    /// Return size of given dimension
    std::size_t size(std::size_t dim) const
    { return PETScBaseMatrix::size(dim); }

    /// Description
    virtual std::string str(bool verbose) const;

  };

}

#endif

#endif
