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
      C_V = std::shared_ptr<void>();
      auto cv = C_V.get();
      cv = CVodeCreate(CV_ADAMS,CV_NEWTON);
    }

    /// Create CVode object with size N
    CVode(MPI_Comm comm, int lmm, int iter)
    {
      C_V = std::shared_ptr<void>();
      auto cv = C_V.get();
      cv = CVodeCreate(lmm, iter);
    }

    /// Destructor 
    ~CVode()
    {
    }

    //--- Implementation of CVode functions

    void init(std::shared_ptr<GenericVector> u0, double t, double atol, double rtol)
    {
//      CVodeInit(C_V,f,t,u0.nvector());
      CVodeSStolerances(C_V.get(),atol, rtol);
    }

    double step(double dt)
    {
    }

    

//    void *CVodeCreate(int lmm, int iter)
//    {
//        void* cv = CVodeCreate(lmm,iter);
//	return cv;	
//   }

    //int SUNDIALSCVodeInit(void* cv,
    //                   int (f)(double,N_Vector,N_Vector, void *),
    //                   double t,
    //                   dolfin::SUNDIALSNVector n0)
    //{
	
//	int flag = CVodeInit(cv,f,t,n0.nvector());
//	return flag;
//    }
    //-----------------------------------------------------------------------------

    /// Assignment operator
    const CVode& operator= (const CVode& x)
    { return *this; }

  private:

    // Pointer to SUNDIALS struct
    std::shared_ptr<void> C_V;

  };

}

#endif

#endif
