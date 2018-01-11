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

#include <dolfin/la/SUNDIALSNVector.h>

#include <cvode/cvode.h>
#include <cvode/cvode_impl.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_iterative.h>

namespace dolfin
{
  /// Wrapper class to SUNDIALS CVODE
  class CVode
  {
  public:

    // These enums are used by PYBIND11 to map the definitions from C
    enum LMM { cv_bdf = CV_BDF, cv_adams = CV_ADAMS };

    enum ITER { cv_functional = CV_FUNCTIONAL, cv_newton = CV_NEWTON };

    /// Constructor
    CVode(LMM cv_lmm, ITER cv_iter);

    /// Destructor
    virtual ~CVode();

    /// Initialise CVode
    void init(std::shared_ptr<GenericVector> u0, double atol, double rtol, long int mxsteps = 0);

    /// Advance time by timestep dt
    double step(double dt);

    /// Get current time
    double get_time() const;

    /// Set the current time
    void set_time(double t0);

    /// Overloaded function for time derivatives of u at time t
    /// Given the vector u, at time t, provide the time derivative udot.
    virtual void derivs(double t, std::shared_ptr<const GenericVector> u,
                        std::shared_ptr<GenericVector> udot);

    /// Given the values (t, y, fy, v), compute Jv = (df/dy)v
    virtual int jacobian(std::shared_ptr<const GenericVector> v,
                         std::shared_ptr<GenericVector> Jv,
                         double t, std::shared_ptr<const GenericVector> y,
                         std::shared_ptr<const GenericVector> fy);

    /// Document
    virtual int jacobian_setup(double t,
                               std::shared_ptr<GenericVector> Jv,
                               std::shared_ptr<GenericVector> y);

    /// Overloaded reconditioner solver function
    virtual int psolve(double tn, std::shared_ptr<const GenericVector>y,
                       std::shared_ptr<const GenericVector> fy,
                       std::shared_ptr<const GenericVector> r,
                       std::shared_ptr<GenericVector> z,
                       double gamma, double delta, int lr);

    /// FIXME: document
    std::map<std::string,double> statistics();

  private:
    // Internal callback from CVode to get time derivatives - passed on to derivs (above)
    static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

    // FIXME: document
    static int f_jac_setup(double t, N_Vector y, N_Vector fy, void *user_data);

    // FIXME: document
    static int f_jac(N_Vector u, N_Vector fu, double t, N_Vector y, N_Vector fy, void* , N_Vector tmp);

    // FIXME: document
    static int prec_solve(double, N_Vector, N_Vector, N_Vector, N_Vector, double, double, int, void*);

    // Vector of values - wrapper around dolfin::GenericVector
    std::shared_ptr<SUNDIALSNVector> _u;

    // SUNDIALS Linear Solver
    SUNLinearSolver _ls;

    // Current time
    double _t;

    // Pointer to CVode memory struct
    void* _cvode_mem;

    // Remember iter method between constructor and init
    ITER _cv_iter;
  };

}

#endif

#endif
