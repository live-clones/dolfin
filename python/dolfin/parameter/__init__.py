
import dolfin.cpp as cpp

#  Extend cpp.Parameters with a __getitem__ method
def __getitem__(self, key):
    if self.has_parameter(key):
        p = self._get_parameter(key)
        return p.value()
    elif self.has_parameter_set(key):
        p = self._get_parameter_set(key)
        np = cpp.parameter.Parameters(p)
        return np
    else:
        raise RuntimeError("invalid parameter")

# Extend the cpp.parameter.Parameters class and clean-up
cpp.parameter.Parameters.__getitem__ = __getitem__
del __getitem__
