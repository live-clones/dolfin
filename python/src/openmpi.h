
#ifndef _DOLFIN_PYBIND11_OPENMPI
#define _DOLFIN_PYBIND11_OPENMPI

// Custom type caster for OpenMPI MPI_Comm
#ifdef OPEN_MPI
#include <mpi.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace pybind11 {
  namespace detail {

    template <> class type_caster<ompi_communicator_t>
    {
    public:
      PYBIND11_TYPE_CASTER(ompi_communicator_t *, _("ompi_communicator_t"));

      // From Python to C++
      bool load(handle src, bool)
      {
        std::uintptr_t v = PyLong_AsUnsignedLong(src.ptr());

        if (PyErr_Occurred()) return false;
        value = reinterpret_cast<ompi_communicator_t *>(v);

        return true;
      }

      // From C++ to Python
      static handle cast(const ompi_communicator_t * const &src, return_value_policy /*policy*/, handle /*parent*/)
      {
        return py::cast(reinterpret_cast<std::uintptr_t>(src));
      }

      operator ompi_communicator_t*()
      {
        return value;
      }
    };
  }
}
// end namespace

#endif
#endif
