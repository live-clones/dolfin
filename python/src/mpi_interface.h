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

#ifndef _DOLFIN_PYBIND11_MPI
#define _DOLFIN_PYBIND11_MPI

#ifdef HAS_MPI
#include <string>
#include <mpi.h>
#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <boost/optional.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//#ifdef HAS_MPI4PY
//#include <mpi4py/mpi4py.h>
//#endif

// Tools for managing MPI communicators

namespace py = pybind11;

namespace pybind11 { namespace detail {
    template <typename T>
    struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};
}}


// `boost::variant` as an example -- can be any `std::variant`-like container
namespace pybind11
{
  namespace detail
  {
    template <typename... Ts>
      struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};

    // Specifies the function used to visit the variant -- `apply_visitor` instead of `visit`
    template <>
      struct visit_helper<boost::variant>
    {
      template <typename... Args>
        static auto call(Args &&...args)
        -> decltype(boost::apply_visitor(std::forward<Args>(args)...))
      {
        return boost::apply_visitor(std::forward<Args>(args)...);
      }
    };
  }
} // namespace pybind11::detail


namespace dolfin_wrappers
{
  // Helper struct for transforming unwrapped DOLFIN MPI communicators
  // in mpi4py communicators
  //struct mpi_communicator
  //{
  //  MPI_Comm comm;
  //};
}

namespace pybind11
{
  namespace detail
  {

//     // Caster used to convert MPI_Comm to a mpi4py communicator
//     template <> class type_caster<dolfin_wrappers::mpi_communicator>
//     {
//     public:
//       PYBIND11_TYPE_CASTER(dolfin_wrappers::mpi_communicator, _("mpi_communicator"));

//       // Helper from pybind to C++
//       bool load(handle src, bool)
//       {
//         //std::cout << "***Here I am" << std::endl;
//         return true;
//       }

//       // From C++ to Python
//       static handle cast(const dolfin_wrappers::mpi_communicator &src, return_value_policy, handle)
//       {
//         //std::cout << "cpp to Py for comm struct" << std::endl;
//         #ifdef HAS_MPI4PY
//         return PyMPIComm_New(src.comm);
//         #else
//         std::cerr << "Cannot convert MPI communicator to a mpi4py communicator. DOLFIN must be enabled \
// with mpi4py for this functionality" << std::endl;
//         #endif
//         return nullptr;
//       }

//       operator dolfin_wrappers::mpi_communicator()
//       {
//         //std::cout << "****mpi comm op" << std::endl;
//         return value;
//       }

//     };


    // Custom type caster for OpenMPI MPI_Comm, in which MPI_Comm is
    // defined as a typedef of ompi_communicator_t*
    #ifdef OPEN_MPI
    template <> class type_caster<ompi_communicator_t>
    {
    public:
      PYBIND11_TYPE_CASTER(MPI_Comm, _("ompi_communicator_t"));

      // Pass communicator from Python to C++
      bool load(handle src, bool)
      {
        PyObject* obj = src.ptr();
        void* v = PyLong_AsVoidPtr(obj);
        value = reinterpret_cast<MPI_Comm>(v);
        if (PyErr_Occurred())
          return false;

        return true;
      }

      // Cast from C++ to Python (cast to pointer)
      static handle cast(MPI_Comm src, py::return_value_policy policy, handle parent)
      {
        // Return a py::handle (rather than a py::object)
        return py::cast(reinterpret_cast<std::uintptr_t>(src)).release();
      }

      operator MPI_Comm()
      { return value; }
    };
//     #else

//     // Custom type caster for MPI_Comm, where MPI_Comm is defined
//     // as a typedef of int, when mpi4py is available
//     #ifdef HAS_MPI4PY
//     template <> class type_caster<MPI_Comm>
//     {
//     public:
//       PYBIND11_TYPE_CASTER(MPI_Comm, _("MPI_Comm"));

//       // From Python (possibly a mpi4py comm) to C++
//       bool load(handle src, bool)
//       {

//         PyObject* obj = src.ptr();

//         if (PyObject_TypeCheck(obj, &PyMPIComm_Type))
//         {
//           MPI_Comm *comm_p = PyMPIComm_Get(obj);
//           value = *comm_p;
//           if (PyErr_Occurred())
//             return false;
//         }
//         else if (PyObject_TypeCheck(obj, &PyLong_Type))
//         {
//           //std::cout << "In caster" << std::endl;
//           value = PyLong_AsLong(obj);
//           if (PyErr_Occurred())
//             return false;
//         }
//         else
//           std::cerr << "MPI communicator (MPI_Comm) type is unknown." << std::endl;

//         return true;
//       }

//       // Cast from C/C++ communicator to Python. Cannot return mpi4py
//       // caster because we cannot distinguish between MPI_Comm and int
//       static handle cast(const MPI_Comm &src, return_value_policy, handle)
//       { return PyLong_FromLong(src); }

//       operator MPI_Comm()
//       { return value; }
//     };
//     #endif
    #endif

  }
}

#endif
#endif
