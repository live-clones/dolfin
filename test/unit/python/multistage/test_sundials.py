#!/usr/bin/python3

from dolfin import *

class MyCVode(CVode):

    # Overload the derivs method
    def derivs(self, t, u, udot):
        udot[:] = -u[:]

def test_sundials():

    phi = Vector(mpi_comm_world(),20)
    phi[:] = 1.0
    cv = MyCVode()
    cv.init(phi, 1e-6, 1e-6)

    nstep = 20
    dt = 0.1
    for i in range(nstep):
        t = cv.step(dt)
        
    assert (exp(-t)-phi.max())<1e-5
