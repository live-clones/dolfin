
import dolfin.cpp as cpp

class Parameters(cpp.parameter.Parameters):
    def __getitem__(self, key):
        if self.has_parameter(key):
            p = self._get_parameter(key)
            return p.value()
        elif self.has_parameter_set(key):
            p = self._get_parameter_set(key)
            np = Parameters(p)
            return np
        else:
            raise RuntimeError("invalid parameter")
