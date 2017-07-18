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

#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/nls/NewtonSolver.h>
#include <dolfin/nls/PETScSNESSolver.h>
#include <dolfin/nls/PETScTAOSolver.h>
#include <dolfin/nls/TAOLinearBoundSolver.h>
#include <dolfin/nls/NonlinearProblem.h>

#include "../mpi_interface.h"

namespace py = pybind11;

namespace dolfin_wrappers
{

  void nls(py::module& m)
  {
    py::class_<dolfin::NewtonSolver>(m, "NewtonSolver")
      .def(py::init<MPI_Comm>());

#ifdef HAS_PETSC
    py::class_<dolfin::PETScSNESSolver>(m, "PETScSNESSolver")
      .def(py::init<MPI_Comm, std::string>());

    py::class_<dolfin::TAOLinearBoundSolver>(m, "TAOLinearBoundSolver")
      .def(py::init<MPI_Comm>());

    py::class_<dolfin::PETScTAOSolver>(m, "PETScTAOSolver")
    .def("__init__",
         [](dolfin::PETScTAOSolver &instance,
            MPI_Comm comm, std::string tao_type,
            std::string ksp_type, std::string pc_type)
         { new (&instance) dolfin::PETScTAOSolver(comm, tao_type, ksp_type, pc_type); },
         py::arg("comm"), py::arg("tao_type")="default",
         py::arg("ksp_type")="default", py::arg("pc_type")="default");
#endif

    // dolfin::NonlinearProblem
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

    py::class_<dolfin::NonlinearProblem, PyNonlinearProblem>(m, "NonlinearProblem")
      .def(py::init<>())
      .def("F", &dolfin::NonlinearProblem::F)
      .def("J", &dolfin::NonlinearProblem::J)
      .def("form", (void (dolfin::NonlinearProblem::*)(dolfin::GenericMatrix&, dolfin::GenericMatrix&,
                                                       dolfin::GenericVector&, const dolfin::GenericVector&))
                    &dolfin::NonlinearProblem::form);

  }

}
