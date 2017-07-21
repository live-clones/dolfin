import dolfin.cpp as cpp

def la_index_dtype():
    """Return the numpy dtype equivalent to the type of la_index"""
    from numpy import intc, int64
    return intc if cpp.common.sizeof_la_index() == 4 else int64


def _unwrap_common_la(x):
    y = x.instance()
    return as_backend_type(y)


def as_backend_type(x):
    converters = [(cpp.la.has_type_vector, cpp.la.as_type_vector),
                  (cpp.la.has_type_matrix, cpp.la.as_type_matrix),
                  (cpp.la.has_type_eigen_vector, cpp.la.as_type_eigen_vector),
                  (cpp.la.has_type_eigen_matrix, cpp.la.as_type_eigen_matrix)]

    converters.extend([(cpp.la.has_type_petsc_vector, cpp.la.as_type_petsc_vector),
                       (cpp.la.has_type_petsc_matrix, cpp.la.as_type_petsc_matrix)])

    for type_check, converter in converters:
        if type_check(x):
            return converter(x)
