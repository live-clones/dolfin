// Copyright (C) 2017 Chris Richardson and Chris Hadjigeorgiou
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

#ifndef __DOLFIN_C_VODE_H
#define __DOLFIN_C_VODE_H

#ifdef HAS_SUNDIALS

#include <cvode/cvode.h>
#include <dolfin/la/SUNDIALSNVector.h>

namespace dolfin
{
  class CVode
  {
  public:

    /// Constructor
    CVode()
    {
      // Create CVode memory block
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      dolfin_assert(cvode_mem);

      // Point user_data back to this object
      // (for use in "f" below)
      int flag = CVodeSetUserData(cvode_mem, (void *)this);
      dolfin_assert(flag == 0);
    }

    /// Destructor
    ~CVode()
    {
      CVodeFree(&cvode_mem);
    }

    /// Initialise CVode
    void init(std::shared_ptr<GenericVector> u0, double atol, double rtol);

    /// Advance time by timestep dt
    double step(double dt);

    /// Get current time
    double time()
    { return t; }

    /// Overloaded function for time derivatives of u at time t
    /// Given the vector u, at time t, provide the time derivative udot.
    virtual void derivs(double t, std::shared_ptr<GenericVector> u,
                        std::shared_ptr<GenericVector> udot);

  private:
    // Internal callback from CVode to get time derivatives - passed on to derivs (above)
    static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

    // Vector of values - wrapper around dolfin::GenericVector
    std::shared_ptr<SUNDIALSNVector> _u;

    // Current time
    double t;

    // Pointer to CVode memory struct
    void *cvode_mem;

  };

}

#endif

#endif
