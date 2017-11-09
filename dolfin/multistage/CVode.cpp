// Copyright (C) 2017 Chris Richardson
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

#ifdef HAS_SUNDIALS

#include <cmath>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include <dolfin/common/MPI.h>
#include <dolfin/la/SUNDIALSNVector.h>
#include <dolfin/log/log.h>

#include <cvode/cvode.h>
#include <cvode/cvode_bandpre.h>
#include <cvode/cvode_impl.h>
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_iterative.h>



#include "CVode.h"

using namespace dolfin;

/// Constructor
CVode::CVode(int cv_lmm, int cv_iter) : t(0.0)
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
CVode::~CVode()
{
  CVodeFree(&cvode_mem);
}

//-----------------------------------------------------------------------------
void CVode::init(std::shared_ptr<GenericVector> u0, double atol, double rtol, long int mxsteps)
{
  dolfin_assert(cvode_mem);

  auto fu = std::shared_ptr<GenericVector>();
  int NEQ = 20*20;

  // Make a sundials n_vector sharing data with u0
  _u = std::make_shared<SUNDIALSNVector>(u0);
  SUNLinearSolver sunls = SUNSPGMR(_u->nvector(), PREC_LEFT,0);

  // Initialise
  std::cout << "Initialising with t = " << t << "\n";

  int flag = CVodeInit(cvode_mem, f, t, _u->nvector());
  dolfin_assert(flag == CV_SUCCESS);

  flag = CVodeSStolerances(cvode_mem, rtol, atol);
  dolfin_assert(flag == CV_SUCCESS);

  CVodeSetMaxNumSteps(cvode_mem, mxsteps);


  if(cv_iter == CV_NEWTON)
  {
    dolfin_debug("Initialising Newtonian solver");
    ls = std::unique_ptr<_generic_SUNLinearSolver>(new _generic_SUNLinearSolver());
    ls->ops = sunls->ops;
    ls->content = sunls->content;
    flag = CVSpilsSetLinearSolver(cvode_mem,ls.get());
    dolfin_assert(flag == CV_SUCCESS);

    /* Call CVBandPreInit to initialize band preconditioner */
    flag = CVSpilsSetPreconditioner(cvode_mem, NULL, CVode::PSolve);
    dolfin_assert(flag == CV_SUCCESS);
    flag = CVSpilsSetJacTimes(cvode_mem, NULL, fJac);
    dolfin_assert(flag == CV_SUCCESS);
  }

}
//-----------------------------------------------------------------------------
double CVode::step(double dt)
{

  double tout = t + dt;
  std::cout << "t_in = " << t << ", dt = " << dt << " " << ((CVodeMem) cvode_mem)->cv_h << std::endl;
  int flag = ::CVode(cvode_mem, tout, _u->nvector(), &t, CV_NORMAL);
  dolfin_assert((flag == CV_SUCCESS) || (flag == CV_TSTOP_RETURN) || (flag == CV_ROOT_RETURN));
  
  std::cout << "t_out = " << t;

  return t;
}
//-----------------------------------------------------------------------------
double CVode::get_time() const
{
  return t;
}
//-----------------------------------------------------------------------------
void CVode::set_time(double t0)
{
  t = t0;
}
//-----------------------------------------------------------------------------
void CVode::derivs(double t, std::shared_ptr<GenericVector> u,
                   std::shared_ptr<GenericVector> udot)
{
  dolfin_error("CVode.cpp",
               "form time derivative",
               "This function should be overloaded");
}
//-----------------------------------------------------------------------------
int CVode::Jacobian(std::shared_ptr<GenericVector> v,
                          std::shared_ptr<GenericVector> Jv,
                          double t, std::shared_ptr<GenericVector> y,
                          std::shared_ptr<GenericVector> fy)
{
  dolfin_error("CVode.cpp",
	       "Jacobian function",
	       "This function should be overloaded");
  return 0;
}
//-----------------------------------------------------------------------------
int CVode::JacobianSetup(double t,
                          std::shared_ptr<GenericVector> Jv,
                          std::shared_ptr<GenericVector> y)
{
  dolfin_error("CVode.cpp",
	       "Jacobian setup function",
	       "This function should be overloaded");
  return 0;
}


//-----------------------------------------------------------------------------
int CVode::fJacSetup(double t, N_Vector y, N_Vector fy, void *user_data)
{

  CVode* cv = static_cast<CVode*>(user_data);

  auto yvec = static_cast<const SUNDIALSNVector*>(y->content)->vec();
  auto fyvec = static_cast<SUNDIALSNVector*>(fy->content)->vec();
  
  cv->JacobianSetup(t, yvec, fyvec);
  return 0;
}
//-----------------------------------------------------------------------------
int CVode::fJac(N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void* user_data, N_Vector tmp)
{

  CVode* cv = static_cast<CVode*>(user_data);

  auto vvec = static_cast<const SUNDIALSNVector*>(v->content)->vec();
  auto Jvvec = static_cast<SUNDIALSNVector*>(Jv->content)->vec();

  auto yvec = static_cast<const SUNDIALSNVector*>(y->content)->vec();
  auto fyvec = static_cast<SUNDIALSNVector*>(fy->content)->vec();
  auto tmpvec = static_cast<SUNDIALSNVector*>(tmp->content)->vec();

  cv->Jacobian(vvec,Jvvec,t,yvec,fyvec);
  return 0;
}
//-----------------------------------------------------------------------------
int CVode::f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  // f is a static function (from C), so need to get pointer to "this" object
  // passed though in user_data
  CVode* cv = static_cast<CVode*>(user_data);

  auto uvec = static_cast<const SUNDIALSNVector*>(u->content)->vec();
  auto udotvec = static_cast<SUNDIALSNVector*>(udot->content)->vec();

  // Callback to actually calculate the derivatives (user function)
  cv->derivs(t, uvec, udotvec);

  return 0;
}

/* Preconditioner solve routine */
/* TODO: Create as virtual function */
int CVode::PSolve(double tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  double gamma, double delta, int lr, void *user_data)
{
  return(0);
}


//-----------------------------------------------------------------------------
std::map<std::string, double> CVode::statistics()
{
  std::map<std::string, double> stats;
  auto cv = static_cast<const CVodeMem>(cvode_mem);

  stats["Steps"] = cv->cv_nst;
  stats["RHSEvals"] = cv->cv_nfe;
  stats["LinSolvSetups"] = cv->cv_nsetups;
  stats["ErrTestFails"] = cv->cv_netf;
  stats["LastOrder"] = cv->cv_qu;
  stats["CurrentOrder"] = cv->cv_next_q;
  stats["StabLimOrderReds"] = cv->cv_nor;
  stats["ActualInitStep"] = cv->cv_h0u;
  stats["LastStep"] = cv->cv_hu;
  stats["CurrentStep"] = cv->cv_next_h;
  stats["CurrentTime"] = cv->cv_tn;
  stats["TolScaleFactor"] = cv->cv_tolsf;
  stats["NumGEvals"] = cv->cv_nge;
  stats["NumNonlinSolvIter"] = cv->cv_nni;
  stats["NumNonlinSolvConvFails"] = cv->cv_ncfn;

  return stats;
}
#endif
