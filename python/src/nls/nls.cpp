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

#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <dolfin/parameter/Parameters.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/PETScObject.h>
#include <dolfin/nls/NewtonSolver.h>
#include <dolfin/nls/PETScSNESSolver.h>
#include <dolfin/nls/PETScTAOSolver.h>
#include <dolfin/nls/TAOLinearBoundSolver.h>
#include <dolfin/nls/NonlinearProblem.h>
#include <dolfin/nls/OptimisationProblem.h>

#include "../mpi_interface.h"

namespace py = pybind11;

namespace dolfin_wrappers
{

  void nls(py::module& m)
  {
    py::class_<dolfin::NewtonSolver, std::shared_ptr<dolfin::NewtonSolver>,
               dolfin::Variable>(m, "NewtonSolver")
      .def(py::init<MPI_Comm>());

#ifdef HAS_PETSC
    py::class_<dolfin::PETScSNESSolver, std::shared_ptr<dolfin::PETScSNESSolver>,
               dolfin::PETScObject>(m, "PETScSNESSolver")
      .def(py::init<MPI_Comm, std::string>());

    py::class_<dolfin::TAOLinearBoundSolver>(m, "TAOLinearBoundSolver")
      .def(py::init<MPI_Comm>())
      .def(py::init<std::string, std::string, std::string>(), py::arg("method")="default",
           py::arg("ksp_type")="default", py::arg("pc_type")="default")
      .def("solve", (std::size_t (dolfin::TAOLinearBoundSolver::*)
                     (const dolfin::GenericMatrix&, dolfin::GenericVector&,
                      const dolfin::GenericVector&, const dolfin::GenericVector&,
                      const dolfin::GenericVector&))
           &dolfin::TAOLinearBoundSolver::solve);

    py::class_<dolfin::PETScTAOSolver, std::shared_ptr<dolfin::PETScTAOSolver>, dolfin::PETScObject>(m, "PETScTAOSolver")
      .def(py::init<>())
      .def("__init__",
           [](dolfin::PETScTAOSolver &instance,
              MPI_Comm comm, std::string tao_type,
              std::string ksp_type, std::string pc_type)
           { new (&instance) dolfin::PETScTAOSolver(comm, tao_type, ksp_type, pc_type); },
           py::arg("comm"), py::arg("tao_type")="default",
           py::arg("ksp_type")="default", py::arg("pc_type")="default")
      .def_readwrite("parameters", &dolfin::PETScTAOSolver::parameters)
        .def("solve", (std::pair<std::size_t, bool> (dolfin::PETScTAOSolver::*)(dolfin::OptimisationProblem&, dolfin::GenericVector&))
             &dolfin::PETScTAOSolver::solve)
        .def("solve", (std::pair<std::size_t, bool> (dolfin::PETScTAOSolver::*)(dolfin::OptimisationProblem&, dolfin::GenericVector&,
                                                                                const dolfin::GenericVector&, const dolfin::GenericVector&))
             &dolfin::PETScTAOSolver::solve);
    //.def("solve", [](dolfin::PETScTAOSolver& self, dolfin::OptimisationProblem& problem, dolfin::GenericVector& x)
    //     { auto val = self.solve(problem, x); return py::make_tuple(val.first, val.second); })
      //ef("solve", [](dolfin::PETScTAOSolver& self, dolfin::OptimisationProblem& problem, dolfin::GenericVector& x,
      //              const dolfin::GenericVector& lb, const dolfin::GenericVector& ub)
      //   {
      //     //std::cout << "Here I am" << std::endl;
      //     std::pair<std::size_t, bool> mval = self.solve(problem, x, lb, ub);
      //    return 1;
      //     //return py::make_tuple(val.first, val.second);
      // });

#endif

    // dolfin::NonlinearProblem 'trampoline'
    class PyNonlinearProblem : public dolfin::NonlinearProblem
    {
      using dolfin::NonlinearProblem::NonlinearProblem;

      void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x) override
      { PYBIND11_OVERLOAD_PURE(void, dolfin::NonlinearProblem, J, A, x); }

      void F(dolfin::GenericVector& b, const dolfin::GenericVector& x) override
      { PYBIND11_OVERLOAD_PURE(void, dolfin::NonlinearProblem, F, b, x); }

      void form(dolfin::GenericMatrix& A, dolfin::GenericMatrix& P,
                dolfin::GenericVector& b, const dolfin::GenericVector& x) override
      { PYBIND11_OVERLOAD(void, dolfin::NonlinearProblem, form, A, P, b, x); }

    };

    // dolfin::NonlinearProblem
    py::class_<dolfin::NonlinearProblem, std::shared_ptr<dolfin::NonlinearProblem>, PyNonlinearProblem>(m, "NonlinearProblem")
      .def(py::init<>())
      .def("F", &dolfin::NonlinearProblem::F)
      .def("J", &dolfin::NonlinearProblem::J)
      .def("form", (void (dolfin::NonlinearProblem::*)(dolfin::GenericMatrix&, dolfin::GenericMatrix&,
                                                       dolfin::GenericVector&, const dolfin::GenericVector&))
                    &dolfin::NonlinearProblem::form);


    // dolfin::OptimizationProblem 'trampoline
    class PyOptimisationProblem : public dolfin::OptimisationProblem
    {
      // See https://github.com/pybind/pybind11/issues/250

      using dolfin::OptimisationProblem::OptimisationProblem;
      //using _binder_base_ = dolfin::OptimisationProblem;
      //using dolfin::OptimisationProblem;

      //double f(const dolfin::GenericVector& x) override
      //{
      //  //PYBIND11_OVERLOAD_PURE(double, dolfin::OptimisationProblem, f, x);
      //  PYBIND11_OVERLOAD_INT(double, dolfin::OptimisationProblem, "f", &x);
      //  return dolfin::OptimisationProblem::f(x);
      //  pybind11::pybind11_fail("Tried to call pure virtual function");
      //}

      double f(const dolfin::GenericVector& x) override
      {
        pybind11::gil_scoped_acquire gil;
        pybind11::function overload = pybind11::get_overload(static_cast<const dolfin::OptimisationProblem *>(this), "f");
        if (overload)
        {
          auto o = overload.operator()<pybind11::return_value_policy::reference>(x);
          if (pybind11::detail::cast_is_temporary_value_reference<double>::value)
          {
            static pybind11::detail::overload_caster_t<double> caster;
            return pybind11::detail::cast_ref<double>(std::move(o), caster);
          }
          else return pybind11::detail::cast_safe<double>(std::move(o));
        }
        pybind11::pybind11_fail("Tried to call pure virtual function \"AAA::pv_foo\"");
      }

      void F(dolfin::GenericVector& b, const dolfin::GenericVector& x) override
      { PYBIND11_OVERLOAD_PURE(void, dolfin::OptimisationProblem, F, b, x); }

      void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x) override
      { PYBIND11_OVERLOAD_PURE(void, dolfin::OptimisationProblem, J, A, x); }

    };

    // dolfin::OptimizationProblem
    py::class_<dolfin::OptimisationProblem, std::shared_ptr<dolfin::OptimisationProblem>,
               PyOptimisationProblem, dolfin::NonlinearProblem>(m, "OptimisationProblem")
      .def(py::init<>())
      .def(py::init<const dolfin::OptimisationProblem&>())
      .def("f", &dolfin::OptimisationProblem::f)
      .def("F", &dolfin::OptimisationProblem::F)
      .def("J", &dolfin::OptimisationProblem::J);

  }

}
