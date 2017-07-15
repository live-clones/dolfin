
import dolfin.cpp as cpp

class Parameters(cpp.parameter.Parameters):
    def __getitem__(self, key):
        p = self._get_parameter(key)
        return p.value()
