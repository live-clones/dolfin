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

#ifndef _DOLFIN_PYBIND11_PETSC
#define _DOLFIN_PYBIND11_PETSC

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef HAS_PETSC4PY
#include <petsc4py/petsc4py.h>

namespace py = pybind11;

namespace pybind11
{
  namespace detail
    {
    template <> class type_caster<_p_KSP>
    {
    public:
      PYBIND11_TYPE_CASTER(KSP, _("ksp"));

      // Pass communicator from Python to C++
      bool load(handle src, bool)
      {
        // FIXME: check reference counting
        std::cout << "Py to c++" << std::endl;
        value = PyPetscKSP_Get(src.ptr());
        /*
        PyObject* obj = src.ptr();
        void* v = PyLong_AsVoidPtr(obj);
        value = reinterpret_cast<MPI_Comm>(v);
        if (PyErr_Occurred())
          return false;
        */

        return true;
      }

      // Cast from C++ to Python (cast to pointer)
      static handle cast(KSP src, py::return_value_policy policy, handle parent)
      {
        // FIXME: check reference counting
        std::cout << "C++ to Python" << std::endl;
        return py::handle(PyPetscKSP_New(src));
      }

      operator KSP()
      { return value; }
    };

    template <> class type_caster<_p_DM>
    {
    public:
      PYBIND11_TYPE_CASTER(DM, _("DM"));

      // Pass communicator from Python to C++
      bool load(handle src, bool)
      {
        // FIXME: check reference counting
        std::cout << "Py to C++ (DM)" << std::endl;
        value = PyPetscDM_Get(src.ptr());
        return true;
      }

      // Cast from C++ to Python (cast to pointer)
      static handle cast(DM src, py::return_value_policy policy, handle parent)
      {
        // FIXME: check reference counting
        std::cout << "C++ to Python (DM)" << std::endl;
        return py::handle(PyPetscDM_New(src));
      }

      operator DM()
      { return value; }
    };
  }
}

#endif
#endif
