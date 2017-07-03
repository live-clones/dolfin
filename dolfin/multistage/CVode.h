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
  /// Wrapper class to SUNDIALS CVODE
  class CVode
  {
  public:

    /// Constructor
    CVode(int cv_lmm = CV_BDF, int cv_iter = CV_ADAMS) : t(0.0)
    {
      
      this->cv_lmm = cv_lmm;
      this->cv_iter = cv_iter;
      // Create CVode memory block
      cvode_mem = CVodeCreate(cv_lmm, cv_iter);
      dolfin_assert(cvode_mem);

      // Point user_data back to this object
      // (for use in "f" below)
      int flag = CVodeSetUserData(cvode_mem, (void *)this);
      dolfin_assert(flag == 0);
    }

    /// Destructor
    virtual ~CVode()
    {
      CVodeFree(&cvode_mem);
    }

    /// Initialise CVode
    void init(std::shared_ptr<GenericVector> u0, double atol, double rtol);

    /// Advance time by timestep dt
    double step(double dt);

    /// Get current time
    double get_time() const
    { return t; }

    /// Set the current time
    void set_time(double t0)
    { t = t0; }

    /// Overloaded function for time derivatives of u at time t
    /// Given the vector u, at time t, provide the time derivative udot.
    virtual void derivs(double t, std::shared_ptr<GenericVector> u,
                        std::shared_ptr<GenericVector> udot);

    /// Overloaded Jabocian function
    virtual int Jacobian(std::shared_ptr<GenericVector> v,
                          std::shared_ptr<GenericVector> Jv,
   		          double t, std::shared_ptr<GenericVector> y,
                          std::shared_ptr<GenericVector> fy);


    std::map<std::string,double> statistics();

  private:
    // Internal callback from CVode to get time derivatives - passed on to derivs (above)
    static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

    static int fJac(N_Vector u, N_Vector fu, double t, N_Vector x, N_Vector y, void* , N_Vector z);

    // Vector of values - wrapper around dolfin::GenericVector
    std::shared_ptr<SUNDIALSNVector> _u;

    // Current time
    double t;
    int cv_lmm, cv_iter;

    // Pointer to CVode memory struct
    void *cvode_mem;

  };

}

#endif

#endif
