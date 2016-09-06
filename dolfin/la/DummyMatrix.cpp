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


#include "DummyMatrix.h"

#include "TensorLayout.h"
#include "EigenFactory.h"

#include <dolfin/common/MPI.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DummyMatrix::DummyMatrix() : DummyMatrix(MPI_COMM_WORLD)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
DummyMatrix::DummyMatrix(MPI_Comm comm) : comm{comm}, sizes{0,0}
{
  // Do nothing
}
//-----------------------------------------------------------------------------
DummyMatrix::~DummyMatrix()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void DummyMatrix::init(const TensorLayout& tensor_layout)
{
  sizes[0] = tensor_layout.size(0);
  sizes[1] = tensor_layout.size(1);
}
//-----------------------------------------------------------------------------
void DummyMatrix::init_vector(GenericVector& z, std::size_t dim) const
{
}
//-----------------------------------------------------------------------------
GenericLinearAlgebraFactory& DummyMatrix::factory() const
{
  return EigenFactory::instance();
}
//-----------------------------------------------------------------------------
