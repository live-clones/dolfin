#!/usr/bin/python3

from dolfin import *

class MyCVode(CVode):

    # Overload the derivs method
    def derivs(self, t, u, udot):
        print("derivs")
        udot[:] = -u[:]
    def jacobian(self, u, udot, t, y, fy):
        udot[:] = u[:]
        return 0
    def psolve(self, t, u, udot, r, z, gamma, x, y):
        z[:] = r[:]
        print("psolve")
        return 0

def test_sundials():

    if has_sundials():

        d = open("data.txt", "w")
        phi = Vector(mpi_comm_world(), 10)
        phi[:] = 1.0
        cv = MyCVode(CVode.cv_bdf, CVode.cv_newton)
        cv.init(phi, 1e-6, 1e-6)

        nstep = 200
        dt = 0.01
        for i in range(nstep):
            t = cv.step(dt)
            print(cv.statistics())
            d.write("%f %f %f\n"%(t, phi[0][0], exp(-t)))
#            assert (exp(-t)-phi[0])<1e-3
        d.close()
