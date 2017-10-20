// Copyright (C) 2017 Chris Richardson and Garth N. Wells
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

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/multistage/CVode.h>
#include <dolfin/multistage/MultiStageScheme.h>
#include <dolfin/multistage/PointIntegralSolver.h>
#include <dolfin/multistage/RKSolver.h>

namespace py = pybind11;

namespace dolfin_wrappers
{

  void multistage(py::module& m)
  {
    // dolfin::MultiStageScheme
    py::class_<dolfin::MultiStageScheme, std::shared_ptr<dolfin::MultiStageScheme>>
      (m, "MultiStageScheme")
      .def(py::init<std::vector<std::vector<std::shared_ptr<const dolfin::Form>>>,
           std::shared_ptr<const dolfin::Form>,
           std::vector<std::shared_ptr<dolfin::Function>>,
           std::shared_ptr<dolfin::Function>,
           std::shared_ptr<dolfin::Constant>,
           std::shared_ptr<dolfin::Constant>,
           std::vector<double>,
           std::vector<int>,
           unsigned int,
           const std::string,
           const std::string,
           std::vector<std::shared_ptr<const dolfin::DirichletBC>>>())
      .def("order", &dolfin::MultiStageScheme::order);

    // dolfin::RKSolver
    py::class_<dolfin::RKSolver, std::shared_ptr<dolfin::RKSolver>>
      (m, "RKSolver")
      .def(py::init<std::shared_ptr<dolfin::MultiStageScheme>>())
      .def("step_interval", &dolfin::RKSolver::step_interval);

    // dolfin::PointIntegralSolver
    py::class_<dolfin::PointIntegralSolver, std::shared_ptr<dolfin::PointIntegralSolver>>
      (m, "PointIntegralSolver")
      .def(py::init<std::shared_ptr<dolfin::MultiStageScheme>>())
      .def("reset_newton_solver", &dolfin::PointIntegralSolver::reset_newton_solver)
      .def("reset_stage_solutions", &dolfin::PointIntegralSolver::reset_stage_solutions)
      .def("step", &dolfin::PointIntegralSolver::step)
      .def("step_interval", &dolfin::PointIntegralSolver::step_interval);

    #ifdef HAS_SUNDIALS
    //CVode trampoline class for allowing overloading of virtual functions
    class PyCVode : public dolfin::CVode
    {
    public:

      using dolfin::CVode::CVode;

      void derivs(double t, std::shared_ptr<dolfin::GenericVector> u,
                  std::shared_ptr<dolfin::GenericVector> udot)
      { PYBIND11_OVERLOAD_NAME(void, dolfin::CVode, "derivs", derivs,
                                t, u, udot);}
      int Jacobian( std::shared_ptr<dolfin::GenericVector> v,
                    std::shared_ptr<dolfin::GenericVector> Jv,
                    double t,
                    std::shared_ptr<dolfin::GenericVector> y,
                    std::shared_ptr<dolfin::GenericVector> fy)
      { PYBIND11_OVERLOAD_NAME(int, dolfin::CVode, "Jacobian", Jacobian,
                                v, Jv, t, y, fy);}
    };

    //dolfin::CVode
    py::class_<dolfin::CVode, PyCVode, std::shared_ptr<dolfin::CVode>>cvode(m,"CVode");
    cvode
      .def(py::init<int, int>())
      .def("init", &dolfin::CVode::init, py::arg("u0"), py::arg("atol"), 
        py::arg("rtol"), py::arg("mxsteps")=0)
      .def("set_time", &dolfin::CVode::set_time, py::arg("t0"))
      .def("get_time", &dolfin::CVode::get_time)
      .def("step", (double (dolfin::CVode::*)(double))
        &dolfin::CVode::step)
      .def("statistics", &dolfin::CVode::statistics)
      .def("derivs", (void (dolfin::CVode::*)(double,
        std::shared_ptr<dolfin::GenericVector>,
        std::shared_ptr<dolfin::GenericVector>))
        &dolfin::CVode::derivs)
      .def("Jacobian", (int (dolfin::CVode::*)(std::shared_ptr<dolfin::GenericVector>,
        std::shared_ptr<dolfin::GenericVector>, double,
        std::shared_ptr<dolfin::GenericVector>,
        std::shared_ptr<dolfin::GenericVector>))
        &dolfin::CVode::derivs);
    py::enum_<dolfin::CVode::LMM>(cvode,"LMM")
      .value("CV_BDF", dolfin::CVode::LMM::cv_bdf)
      .value("CV_ADAMS", dolfin::CVode::LMM::cv_adams); 
    py::enum_<dolfin::CVode::ITER>(cvode,"ITER")
      .value("CV_FUNCTIONAL", dolfin::CVode::ITER::cv_functional)
      .value("CV_NEWTON", dolfin::CVode::ITER::cv_newton);
    #endif

  }
}
