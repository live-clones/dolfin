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
#include <dolfin/common/constants.h>
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
    m.def("git_commit_hash", &dolfin::git_commit_hash);
    m.def("sizeof_la_index", &dolfin::sizeof_la_index);

    m.attr("DOLFIN_EPS") = DOLFIN_EPS;
    m.attr("DOLFIN_PI") = DOLFIN_PI;
  }

  // Return a MPI communicator from a mpi4py comm (untested)
  /*
  static MPI_Comm mpi4py_to_comm(py::object obj)
  {
    MPI_Comm *comm_p = PyMPIComm_Get(obj.ptr());
    if (!comm_p)
    {
      std::cout << "Cast failed" << std::endl;
      //throw_error_already_set();
    }

    return *comm_p;
  }
  */

  void mpi(py::module& m)
  {
    #ifdef HAS_MPI4PY
    import_mpi4py();
    #endif

    py::class_<dolfin::MPI>(m, "MPI", "MPI utilities")
      #ifdef OPEN_MPI
      .def_property_readonly_static("comm_world", [](py::object)
                                    { return reinterpret_cast<std::uintptr_t>(MPI_COMM_WORLD); })
      .def_property_readonly_static("comm_self", [](py::object)
                                    { return reinterpret_cast<std::uintptr_t>(MPI_COMM_SELF); })
      .def_property_readonly_static("comm_null", [](py::object)
                                   { return reinterpret_cast<std::uintptr_t>(MPI_COMM_NULL); })
       #else
      .def_property_readonly_static("comm_world", [](py::object) { return MPI_COMM_WORLD; })
      .def_property_readonly_static("comm_self", [](py::object) { return MPI_COMM_SELF; })
      .def_property_readonly_static("comm_null", [](py::object) { return MPI_COMM_NULL; })
       #endif
      .def_static("init", [](){ dolfin::SubSystemsManager::init_mpi();})
      .def_static("barrier", &dolfin::MPI::barrier)
      .def_static("rank", &dolfin::MPI::rank)
      .def_static("size", &dolfin::MPI::size)
      .def_static("max", &dolfin::MPI::max<double>)
      .def_static("min", &dolfin::MPI::min<double>)
      .def_static("sum", &dolfin::MPI::sum<double>) ;
      //.def_static("to_comm", &mpi4py_to_comm);
  }


}
