import dolfin.cpp as cpp


def _add(self, incr):
    self._increment(incr)


def _set(self, value):
    self._assign(value)


def __iadd__(self, other):
    """ Add value to Progress """
    self._add(other)
    return self


cpp.log.Progress._add = _add
cpp.log.Progress.__iadd__ = __iadd__
del __iadd__, _add
