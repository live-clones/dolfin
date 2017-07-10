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

#ifndef _DOLFIN_PYBIND11_OPENMPI
#define _DOLFIN_PYBIND11_OPENMPI

// Custom type caster for OpenMPI MPI_Comm, where MPI_Comm is defined
// as a typedef of ompi_communicator_t*

#ifdef OPEN_MPI
#include <mpi.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace pybind11
{
  namespace detail
  {

    template <> class type_caster<ompi_communicator_t>
    {
    public:
      PYBIND11_TYPE_CASTER(MPI_Comm, _("ompi_communicator_t"));

      // From Python to C++
      bool load(handle src, bool)
      {
        std::cout << "Py to C++ (ptr): " << src.ptr() << std::endl;
        std::cout << "Py to C++ (int): " << reinterpret_cast<std::uintptr_t>(src.ptr()) << std::endl;
        void* v = PyLong_AsVoidPtr(src.ptr());

        if (PyErr_Occurred())
          return false;
        value = reinterpret_cast<MPI_Comm>(v);

        return true;
      }

      // From C++ to Python
      static handle cast(const MPI_Comm &src,
                         return_value_policy, handle)
      {
        std::cout << "C++ to Py: " << &(*src) << std::endl;
        std::cout << "C++ to Py: " << reinterpret_cast<std::uintptr_t>(src) << std::endl;
        return py::cast(reinterpret_cast<std::uintptr_t>(src)); }

      operator MPI_Comm()
      { return value; }
    };
  }
}

#endif
#endif
