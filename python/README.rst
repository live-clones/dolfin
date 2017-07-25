Experimental Python wrapping with pybind11
==========================================

In the directory `python/`

1. Get pybind11 development version::

   git clone https://github.com/pybind/pybind11.git
   cd pybind11
   git fetch origin pull/960/head:registered-type-fixes
   git checkout registered-type-fixes

2. Build bindings::

   python3 setup.py build

3. Add to `PYTHONPATH`, e.g.::

   export PYTHONPATH=build/lib.linux-x86_64-3.5/

4. Run test program `poisson.py` in `tests`
