// Copyright (C) 2007 Garth N. Wells
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
// Modified by Chris Hadjigeorgiou: 2017.
//
// First added:  2017-01-31
// Last changed: 2011-02-07

#ifndef __DOLFIN_C_VODE_H
#define __DOLFIN_C_VODE_H

#ifdef HAS_SUNDIALS

#include <string>
#include <utility>
#include <memory>
#include <dolfin/common/types.h>
#include <cvode/cvode.h>
#include <dolfin/la/SUNDIALSNVector.h>

namespace dolfin
{

  template<typename T> class Array;

  class CVode
  {
  public:

    /// Create empty CVode object 
    CVode(MPI_Comm comm=MPI_COMM_WORLD)
    {
      C_V = CVodeCreate(CV_BDF,CV_NEWTON);
    }

    /// Create CVode object with size N
    CVode(MPI_Comm comm, int lmm, int iter)
    {
      C_V = CVodeCreate(lmm, iter);
    }

    /// Destructor 
    ~CVode()
    {
      CVodeFree(&C_V);
    }

    //--- Implementation of CVode functions

    void init(std::shared_ptr<dolfin::SUNDIALSNVector> u0, CVRhsFn f, double t, double atol, double rtol)
    {
      CVodeInit(C_V,f,t,u0->nvector());
      CVodeSStolerances(C_V,atol, rtol);
    }

    double step(std::shared_ptr<dolfin::SUNDIALSNVector> u0, double dt, double *t)
    {
      ::CVode(C_V, dt, u0->nvector(), t, CV_NORMAL);
    }

    void * cv_mem()
    {  return C_V; }
    //-----------------------------------------------------------------------------

    /// Assignment operator
    const CVode& operator= (const CVode& x)
    { return *this; }
    

  private:

    // Pointer to CVode memory struct
    void * C_V;

  };

}

#endif

#endif
