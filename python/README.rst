Experimental Python wrapping with pybind11
==========================================

Using CMake
-----------

1. cmake .
2. make

To set the Python version::

  cmake -DPYTHON_EXECUTABLE:FILEPATH=<path-to-python-executable> .


Using setup.py
--------------

1. Build using::

     python3 setup.py build

   This will build the module in a 'build' directory.

2. Add the build directory to your `PYTHONPATH`.

3. Run the test program `test.py`

The setup script is still rather primitive and will not auto-detect a
number of thing. It assume (for now) that DOLFIN is built with MPI
enabled.
