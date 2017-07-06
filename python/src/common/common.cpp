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
#include <string>
#include <pybind11/pybind11.h>

#include <dolfin/common/MPI.h>
#include <dolfin/common/defines.h>
#include <dolfin/common/SubSystemsManager.h>
#include <dolfin/common/Variable.h>

#include "../openmpi.h"

namespace py = pybind11;

namespace dolfin_wrappers
{
  void common(py::module& m)
  {
    // Variable
    py::class_<dolfin::Variable, std::shared_ptr<dolfin::Variable>>
      (m, "Variable", "Variable base class")
      .def("id", &dolfin::Variable::id)
      .def("name", &dolfin::Variable::name)
      .def("rename", &dolfin::Variable::rename);

    // From dolfin/common/defines.h
    m.def("has_debug", &dolfin::has_debug);
    m.def("has_hdf5", &dolfin::has_hdf5);
    m.def("has_hdf5_parallel", &dolfin::has_hdf5_parallel);
    m.def("has_mpi", &dolfin::has_mpi);
    m.def("has_petsc", &dolfin::has_petsc);
    m.def("has_slepc", &dolfin::has_slepc);
    m.def("git_commit_has", &dolfin::git_commit_hash);

  }

  void mpi(py::module& m)
  {
    // MPI
    // m.attr("comm_world") = MPI_COMM_WORLD;
    // m.attr("comm_self") = MPI_COMM_SELF;

    m.def("init", [](){ dolfin::SubSystemsManager::init_mpi();});
    m.def("my_init", [](){
        int argc = 0;
        std::string s("");
        char* c = const_cast<char *>(s.c_str());
        char **argv = &c;
        MPI_Init(&argc, &argv);
          });

    m.def("comm_world", []() { return MPI_COMM_WORLD; });
    m.def("comm_self", []() { return MPI_COMM_SELF; });

    m.def("rank", &dolfin::MPI::rank);
    m.def("size", &dolfin::MPI::size);
    m.def("max", &dolfin::MPI::max<double>);
    m.def("min", &dolfin::MPI::min<double>);
    m.def("sum", &dolfin::MPI::sum<double>);

  }


}
