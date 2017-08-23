Experimental Python wrapping with pybind11
==========================================

1. Install the development version of pybdind11. Use CMake to install.

   git clone https://github.com/pybind/pybind11.git
   cd pybind11
   mkdir build-dir
   cmake -DPYBIND11_TEST=off -DCMAKE_INSTALL_PREFIX=/path/to/pybind11/install ../
   make install

2. Build bindings::

   export PYBIND11_DIR=/path/to/pybind11/install
   pip3 -v install --user

2. Run test program `poisson.py` in `tests`
