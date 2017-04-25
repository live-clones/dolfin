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

#include <dolfin/common/MPI.h>
#include <dolfin/la/SUNDIALSNVector.h>
#include <dolfin/log/log.h>

#include <vector>
#include <iostream>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

#include "CVode.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void CVode::init(std::shared_ptr<GenericVector> u0, double atol, double rtol)
{
  dolfin_assert(cvode_mem);

  // Make a sundials n_vector sharing data with u0
  _u = std::make_shared<SUNDIALSNVector>(u0);

  // Initialise
  std::cout << "Initialising with t = " << t << "\n";

  int flag = CVodeInit(cvode_mem, f, t, _u->nvector());
  dolfin_assert(flag == CV_SUCCESS);

  flag = CVodeSStolerances(cvode_mem, rtol, atol);
  dolfin_assert(flag == CV_SUCCESS);
}
//-----------------------------------------------------------------------------
double CVode::step(double dt)
{
  //  std::cout << "t_in = " << t;

  double tout = t + dt;
  int flag = ::CVode(cvode_mem, tout, _u->nvector(), &t, CV_NORMAL);
  dolfin_assert(flag == CV_SUCCESS);

  //  std::cout << "t_out = " << t;

  return t;
}
//-----------------------------------------------------------------------------
int CVode::f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  // f is a static function (from C), so need to get pointer to "this" object
  // passed though in user_data
  CVode* cv = static_cast<CVode*>(user_data);

  auto uvec = static_cast<SUNDIALSNVector*>(u->content)->vec();
  auto udotvec = static_cast<SUNDIALSNVector*>(udot->content)->vec();

  // Callback to actually calculate the derivatives (user function)
  cv->derivs(t, uvec, udotvec);

  return 0;
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

#endif
