import dolfin.cpp as cpp

def la_index_dtype():
    """Return the numpy dtype equivalent to the type of la_index"""
    from numpy import intc, int64
    return intc if cpp.common.sizeof_la_index() == 4 else int64
