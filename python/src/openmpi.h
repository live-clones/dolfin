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

#ifdef HAS_MPI
#include <mpi.h>
#include <pybind11/pybind11.h>
#include <Python.h>

#ifdef HAS_MPI4PY
#include <mpi4py/mpi4py.h>
#endif

namespace py = pybind11;

namespace pybind11
{
  namespace detail
  {

    #ifdef OPEN_MPI
    template <> class type_caster<ompi_communicator_t>
    {
    public:
      PYBIND11_TYPE_CASTER(MPI_Comm, _("c_mpi_communicator_t"));

      // From Python to C++
      bool load(handle src, bool)
      {
        PyObject* obj = src.ptr();
        #ifdef HAS_MPI4PY
        if (PyObject_TypeCheck(obj, &PyMPIComm_Type))
        {
          MPI_Comm *comm_p = PyMPIComm_Get(obj);
          value = *comm_p;
          if (PyErr_Occurred())
            return false;
        }
        else
        {
          void* v = PyLong_AsVoidPtr(obj);
          value = reinterpret_cast<MPI_Comm>(v);
          if (PyErr_Occurred())
            return false;
        }
        #else
        void* v = PyLong_AsVoidPtr(obj);
        value = reinterpret_cast<MPI_Comm>(v);
        if (PyErr_Occurred())
          return false;
        #endif

        return true;
      }

      // From C++ to Python
      static handle cast(const MPI_Comm &src,
                         return_value_policy, handle)
      { return py::cast(reinterpret_cast<std::uintptr_t>(src)); }

      operator MPI_Comm()
      { return value; }
    };
    #else
    template <> class type_caster<MPI_Comm>
    {
    public:
      PYBIND11_TYPE_CASTER(MPI_Comm, _("MPI_Comm"));

      // From Python (possibly a mpi4py comm) to C++
      bool load(handle src, bool)
      {
        //std::cout << "Start cast" << std::endl;

        PyObject* obj = src.ptr();
        #ifdef HAS_MPI4PY
        if (PyObject_TypeCheck(obj, &PyMPIComm_Type))
        {
          MPI_Comm *comm_p = PyMPIComm_Get(obj);
          value = *comm_p;
          if (PyErr_Occurred())
            return false;
        }
        else if (PyObject_TypeCheck(obj, &PyLong_Type))
        {
          //std::cout << "In caster" << std::endl;
          value = PyLong_AsLong(obj);
          if (PyErr_Occurred())
            return false;
        }
        else
          std::cerr << "MPI communicator (MPI_Comm) type is unknown." << std::endl;
        #else
        value = PyLong_AsLong(obj);
        if (PyErr_Occurred())
          return false;
        #endif

        return true;
      }

      // Cast from C/C++ communicator to Python
      static handle cast(const MPI_Comm &src, return_value_policy, handle)
      {
        //#if HAS_MPI4PY
        //return PyMPIComm_New(src);
        //return PyMPIComm_New(src);
        //#else
        return PyLong_FromLong(src);
        //#endif
      }

      operator MPI_Comm()
      { return value; }
    };
    #endif

  }
}

#endif
#endif
