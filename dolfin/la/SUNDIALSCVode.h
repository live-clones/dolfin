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

#ifndef __DOLFIN_SUNDIALS_C_VODE_H
#define __DOLFIN_SUNDIALS_C_VODE_H

#ifdef HAS_SUNDIALS

#include <string>
#include <utility>
#include <memory>
#include <dolfin/common/types.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode.h>
#include "DefaultFactory.h"
#include "GenericVector.h"
#include "Vector.h"
#include "SUNDIALSNVector.h"

namespace dolfin
{

  template<typename T> class Array;

  class SUNDIALSCVode
  {
  public:

    /// Create empty CVode object 
    SUNDIALSCVode(MPI_Comm comm=MPI_COMM_WORLD)
    {
    }

    /// Create CVode object with size N
    SUNDIALSCVode(MPI_Comm comm, std::size_t N)
    {
    }

    /// Copy constructor
//    SUNDIALSCVode(const SUNDIALSCVode& x) : vector(x.vector->copy()) {}

    /// Create an SUNDIALSCVode from a GenericVector
//    SUNDIALSCVode(const GenericVector& x) : vector(x.copy())
//    {
//      C_V = std::make_shared<CVodeMem>();
//    }

    //--- Implementation of CVode functions


    void *SUNDIALSCVodeCreate(int lmm, int iter)
    {
        void* cv = CVodeCreate(lmm,iter);
	return cv;	
    }

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
    const SUNDIALSCVode& operator= (const SUNDIALSCVode& x)
    { return *this; }

  private:


    // Pointer to SUNDIALS struct
    std::shared_ptr<void> C_V;


/* Structure containing function pointers to vector operations  */
    struct _generic_N_Vector_Ops ops = {N_VGetVectorID,    //   N_Vector_ID (*N_VGetVectorID)(SUNDIALSCVode);
                                        NULL,    //   CVode    (*N_VClone)(CVode);
                                        NULL,    //   CVode    (*N_VCloneEmpty)(CVode);
                                        NULL,    //   void        (*N_VDestroy)(CVode);
                                        NULL,    //   void        (*N_VSpace)(CVode, long int *, long int *);
                                        NULL,    //   realtype*   (*N_VGetArrayPointer)(CVode);
                                        NULL,    //   void        (*N_VSetArrayPointer)(realtype *, CVode);
                                        NULL,    //   void        (*N_VLinearSum)(realtype, CVode, realtype, CVode, CVode);
                                        N_VConst,          //   void        (*N_VConst)(realtype, CVode);
                                        NULL,    //   void        (*N_VProd)(CVode, CVode, CVode);
                                        NULL,    //   void        (*N_VDiv)(CVode, CVode, CVode);
                                        NULL,    //   void        (*N_VScale)(realtype, CVode, CVode);
                                        NULL,    //   void        (*N_VAbs)(CVode, CVode);
                                        NULL,    //   void        (*N_VInv)(CVode, CVode);
                                        NULL,    //   void        (*N_VAddConst)(CVode, realtype, CVode);
                                        NULL,    //   realtype    (*N_VDotProd)(CVode, CVode);
                                        NULL,    //   realtype    (*N_VMaxNorm)(CVode);
                                        NULL,    //   realtype    (*N_VWrmsNorm)(CVode, CVode);
                                        NULL,    //   realtype    (*N_VWrmsNormMask)(CVode, CVode, CVode);
                                        NULL,    //   realtype    (*N_VMin)(CVode);
                                        NULL,    //   realtype    (*N_VWl2Norm)(CVode, CVode);
                                        NULL,    //   realtype    (*N_VL1Norm)(CVode);
                                        NULL,    //   void        (*N_VCompare)(realtype, CVode, CVode);
                                        NULL,    //   booleantype (*N_VInvtest)(CVode, CVode);
                                        NULL,    //   booleantype (*N_VConstrMask)(CVode, CVode, CVode);
                                        NULL};    //   realtype    (*N_VMinQuotient)(CVode, CVode);
  };


}

#endif

#endif
