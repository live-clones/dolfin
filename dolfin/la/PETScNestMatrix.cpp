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

#ifdef HAS_PETSC

#include "GenericMatrix.h"
#include "GenericVector.h"
#include "PETScMatrix.h"
#include "PETScVector.h"
#include "PETScNestMatrix.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
PETScNestMatrix::PETScNestMatrix()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScNestMatrix::PETScNestMatrix
(std::vector<std::shared_ptr<const GenericMatrix>> mats)
{
  if (mats.size() != 4)
  {
    dolfin_error("PETScNestMatrix.cpp",
                 "create PETScNestmatrix",
                 "Only support 2x2 so far");
  }

  std::vector<Mat> petsc_mats(mats.size());
  MPI_Comm mpi_comm = MPI_COMM_NULL;
  for (std::size_t i = 0; i != mats.size(); ++i)
  {
    std::cout << i << ":";

    if (mats[i])
    {
      petsc_mats[i] = as_type<const PETScMatrix>(mats[i])->mat();
      std::cout << mats[i]->size(0)<< " " << mats[i]->size(1) <<" ";
      // If Mat has been initialised, get mpi_comm
      if (petsc_mats[i])
      {
        if (mpi_comm == MPI_COMM_NULL)
          mpi_comm = mats[i]->mpi_comm();
        else if (mpi_comm != mats[i]->mpi_comm())
        {
          dolfin_error("PETScNestMatrix.cpp",
                       "construct MatNest",
                       "Constituent matrices have different communicators");
        }
      }
    }
    else
      petsc_mats[i] = NULL;
    std::cout << "\n";
  }

  if (mpi_comm == MPI_COMM_NULL)
  {
    dolfin_error("PETScNestMatrix.cpp",
                 "construct MatNest",
                 "All matrices appear to be NULL");
  }

  MatCreateNest(mpi_comm, 2, NULL, 2, NULL, petsc_mats.data(), &_matA);
}
//-----------------------------------------------------------------------------
PETScNestMatrix::~PETScNestMatrix()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::string PETScNestMatrix::str(bool verbose) const
{
  if (verbose)
  {
    PetscViewerSetFormat((PetscViewer)PETSC_VIEWER_DEFAULT,
                         (PetscViewerFormat)PETSC_VIEWER_ASCII_INFO);

    MatView(_matA, (PetscViewer)PETSC_VIEWER_DEFAULT);
  }
  return std::string("PETScNestMatrix");
}
//-----------------------------------------------------------------------------
void PETScNestMatrix::mult(const GenericVector& x, GenericVector& y) const
{
  dolfin_assert(_matA);

  const PETScVector& xx = as_type<const PETScVector>(x);
  PETScVector& yy = as_type<PETScVector>(y);

  if (size(1) != xx.size())
  {
    dolfin_error("PETScNestMatrix.cpp",
                 "compute matrix-vector product with PETSc matrix",
                 "Non-matching dimensions for matrix-vector product");
  }

  if (size(0) != yy.size())
  {
    dolfin_error("PETScNestMatrix.cpp",
                 "compute matrix-vector product with PETSc matrix",
                 "Vector for matrix-vector result has wrong size");
  }

  PetscErrorCode ierr = MatMult(_matA, xx.vec(), yy.vec());
  if (ierr != 0) petsc_error(ierr, __FILE__, "MatMult");
}
//-----------------------------------------------------------------------------


#endif
