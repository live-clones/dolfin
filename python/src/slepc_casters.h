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

#ifndef _DOLFIN_PYBIND11_SLEPC
#define _DOLFIN_PYBIND11_SLEPC

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef HAS_SLEPC
#include <slepceps.h>
// pybind11 casters for SLEPc/slepc4py objects
#ifdef HAS_PYBIND11_SLEPC4PY
#include <slepc4py/slepc4py.h>

// Import slepc4py on demand
#define VERIFY_SLEPC4PY(func)   \
  if (!func)                    \
  {                             \
    if (import_slepc4py() != 0) \
    {                           \
      std::cout << "ERROR: could not import slepc4py!" << std::endl; \
      throw std::runtime_error("Error when importing slepc4py");     \
    }                           \
  }

// Macro for casting between dolfin and petsc4py objects
#define SLEPC_CASTER_MACRO(TYPE, NAME)          \
  template <> class type_caster<_p_##TYPE>      \
    {                                           \
    public:                                     \
      PYBIND11_TYPE_CASTER(TYPE, _(#NAME));     \
      bool load(handle src, bool)               \
      {                                         \
        VERIFY_SLEPC4PY(PySlepc##TYPE##_Get);   \
        if (PyObject_TypeCheck(src.ptr(), &PySlepc##TYPE##_Type) == 0)  \
          return false;                                                 \
        value = PySlepc##TYPE##_Get(src.ptr());                         \
        return true;                                                    \
      }                                                                 \
                                                                        \
      static handle cast(TYPE src, pybind11::return_value_policy policy, handle parent) \
      {                                                                 \
        VERIFY_SLEPC4PY(PySlepc##TYPE##_New);                           \
        return pybind11::handle(PySlepc##TYPE##_New(src));              \
      }                                                                 \
                                                                        \
      operator TYPE()                                                   \
      { return value; }                                                 \
    }
#else
#define SLEPC_CASTER_MACRO(TYPE, NAME)          \
  template <> class type_caster<_p_##TYPE>      \
    {                                           \
    public:                                     \
      PYBIND11_TYPE_CASTER(TYPE, _(#NAME));     \
      bool load(handle src, bool)               \
      {                                         \
        throw std::runtime_error("DOLFIN has not been configured with slepc4py. Accessing underlying SLEPc object requires slepc4py"); \
        return false;                                                   \
      }                                                                 \
                                                                        \
      static handle cast(TYPE src, pybind11::return_value_policy policy, handle parent) \
      {                                                                 \
        throw std::runtime_error("DOLFIN has not been configured with slepc4py. Accessing underlying SLEPc object requires slepc4py"); \
        return handle();                                                \
      }                                                                 \
                                                                        \
      operator TYPE()                                                   \
      { return value; }                                                 \
    }
#endif


namespace pybind11
{
  namespace detail
  {
    SLEPC_CASTER_MACRO(EPS, eps);
  }
}

#undef SLEPC_CASTER_MACRO

#endif
#endif
