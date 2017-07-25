// Copyright (C) 2017 Nate Sime
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
// First added:  2017-07-25


#ifndef __DOLFIN_CONTACTASSEMBLER_H
#define __DOLFIN_CONTACTASSEMBLER_H

#include "AssemblerBase.h"
#include "Assembler.h"

namespace dolfin
{
  class ContactAssembler : Assembler
  {

  public:

    ContactAssembler(const GeometricContact& gc);

    /// Assemble tensor from given form
    ///
    /// @param     A (_GenericTensor_)
    ///         The tensor to assemble.
    /// @param     a (_Form_)
    ///         The form to assemble the tensor from.
    void assemble(GenericTensor& A, const Form& a);

  private:
    void _init_global_tensor(GenericTensor& A, const Form& a);

    const GeometricContact& _gc;

  };
}


#endif
