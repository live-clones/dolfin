#!/usr/bin/python3

from dolfin import *
from dolfin_utils.test import skip_if_not_SUNDIALS

class MyCVode(CVode):

    def derivs(self, t, u, udot):
        udot[:] = -u[:]

    def jacobian(self, v, Jv, t, y, fy):
        Jv[:] = v[:]
        return 0

    def psolve(self, t, u, udot, r, z, gamma, x, y):
        z[:] = r[:]
        return 0

@skip_if_not_SUNDIALS
def test_sundials_adams():

    phi = Vector(MPI.comm_world, 10)
    phi[:] = 1.0

    cv = MyCVode(CVode.LMM.CV_ADAMS, CVode.ITER.CV_FUNCTIONAL)
    cv.init(phi, 1e-7, 1e-7)

    nstep = 200
    dt = 0.01
    for i in range(nstep):
        t = cv.step(dt)
        assert (exp(-t) - phi[0]) < 1e-6

@skip_if_not_SUNDIALS
def test_sundials_newton():

    phi = Vector(MPI.comm_world, 10)
    phi[:] = 1.0

    cv = MyCVode(CVode.LMM.CV_BDF, CVode.ITER.CV_NEWTON)
    cv.init(phi, 1e-7, 1e-7)

    nstep = 200
    dt = 0.01
    for i in range(nstep):
        t = cv.step(dt)
        assert (exp(-t) - phi[0]) < 1e-6

def test_sundials_diffusion():

    mesh = UnitIntervalMesh(MPI.comm_world, 100)
    gaussian = Expression("exp(-pow((x[0]-0.5)/w, 2))", w=0.1, degree=1)
    Q = FunctionSpace(mesh, "CG", 1)
    phi = interpolate(gaussian, Q)

    # Calculate laplacian of phi
    u = TestFunction(Q)
    v = TrialFunction(Q)
    a = u*v*dx
    L = -dot(grad(phi), grad(v))*dx
    A = assemble(a)
    b = assemble(L)
    solver = LUSolver(A)

    class tmp_test(CVode):
        def derivs(self, t, phi, phidot):
            b = assemble(L)
            solver.solve(phidot, b)
        def jacobian(self, v, Jv, t, y, fy):
            Jv[:] = v[:]
            return 0
        def psolve(self, t, u, udot, r, z, gamma, x, y):
            z[:] = r[:]
            return 0

    cv = tmp_test(CVode.LMM.CV_BDF, CVode.ITER.CV_NEWTON)
    cv.init(phi.vector(), 1e-7, 1e-7)

    nstep = 200
    dt = 0.0001
    for i in range(nstep):
        t = cv.step(dt)
        print(t, phi.vector().max())
