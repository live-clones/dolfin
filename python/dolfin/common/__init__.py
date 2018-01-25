import dolfin.cpp as cpp


def to_dict(self):
    """Convert the Parameters to a dict"""
    ret = {}
    for key, value in self.items():
        if isinstance(value, cpp.parameter.Parameters):
            ret[key] = value.to_dict()
        else:
            ret[key] = value
    return ret


cpp.parameter.Parameters.to_dict = to_dict
del to_dict


def _add(self, incr):
    for j in range(incr):
        self._increment(1)


def _set(self, value):
    self._assign(value)


def __iadd__(self, other):
    """ Add value to Progress """
    if isinstance(other, int):
        self._add(other)
    elif isinstance(other, float):
        self._set(other)
    return self


cpp.log.Progress._add = _add
cpp.log.Progress.__iadd__ = __iadd__
del __iadd__, _add
