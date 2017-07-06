
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
      PYBIND11_TYPE_CASTER(MPI_Comm, _("ompi_communicator_t"));

      // From Python to C++
      bool load(handle src, bool)
      {
        void* v = PyLong_AsVoidPtr(src.ptr());

        if (PyErr_Occurred()) return false;
        value = reinterpret_cast<MPI_Comm>(v);

        return true;
      }

      // From C++ to Python
      static handle cast(const MPI_Comm &src,
                         return_value_policy /*policy*/, handle /*parent*/)
      {
        return py::cast(reinterpret_cast<std::uintptr_t>(src));
      }

      operator MPI_Comm()
      {
        return value;
      }
    };
  }
}
// end namespace

#endif
#endif
