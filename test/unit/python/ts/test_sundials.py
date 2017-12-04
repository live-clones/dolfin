#!/usr/bin/python3

from dolfin import *

class MyCVode(CVode):

    # Overload the derivs method
    def derivs(self, t, u, udot):
        udot[:] = -u[:]
    def Jacobian(self, u, udot, t, y, fy):
        return 0
    def psolve(self, t, u, udot, r, z, gamma, x, y):
        return 0

def test_sundials():

    if has_sundials():
      phi = Vector(mpi_comm_world(),20)
      phi[:] = 1.0
      cv = MyCVode(CVode.LMM.CV_BDF,CVode.ITER.CV_NEWTON)
      cv.init(phi, 1e-6, 1e-6)

      nstep = 20
      dt = 0.001
      for i in range(nstep):
        t = cv.step(dt)       
        assert (exp(-t)-phi[i])<1e-3
