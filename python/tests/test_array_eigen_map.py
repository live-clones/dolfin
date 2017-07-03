
import timeit

setup="""
import dolfin_test.cpp as cpp
from dolfin_test.cpp.mesh import SubDomain
w = SubDomain()
"""

t = timeit.timeit("w.test(True, 100000)", setup=setup, number=5000)
print(t)

t = timeit.timeit("w.test(False, 100000)", setup=setup, number=5000)
print(t)
