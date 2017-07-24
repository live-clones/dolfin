import dolfin.cpp as cpp

def la_index_dtype():
    """Return the numpy dtype equivalent to the type of la_index"""
    from numpy import intc, int64
    return intc if cpp.common.sizeof_la_index() == 4 else int64


def as_backend_type(x):
    """Return Matrix and Vector backend instance. Not required for other
    types as pybind11 automatically downcasts objects to the derived
    type.

    """
    if isinstance(x, cpp.la.Vector) or isinstance(x, cpp.la.Matrix):
        return x.instance()
    else:
        return x
