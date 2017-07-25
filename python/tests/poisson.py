from dolfin import *

if cpp.common.has_petsc():
    parameters['linear_algebra_backend'] = 'PETSc'

# Create mesh and refine
mesh = UnitSquareMesh(12, 12)
mesh = refine(mesh)

# Create function space
V = FunctionSpace(mesh, "Lagrange", 1)

# Create a function
w = Function(V)

#xdmf = XDMFFile("a.xdmf")
#xdmf.write(mesh,  XDMFFile.Encoding.ASCII)
#xdmf.write(w, XDMFFile.Encoding.ASCII)

#class Boundary(SubDomain):
#    def inside(self, x, on_boundary):
#        result = (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS)
#        return bool(result)
# boundary = CompiledSubDomain("x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS")

u0 = Constant(0.0)
bc = DirichletBC(V, u0, "x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS")

u = TrialFunction(V)
v = TestFunction(V)
f = CompiledExpression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = CompiledExpression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

assembler = Assembler()

A = Matrix()
assembler.assemble(A, Form(a))

b = Vector()
myform = Form(L)
assembler.assemble(b, myform)

bc.apply(b)
bc.apply(A)

solver = KrylovSolver(A)
solver.solve(w.vector(), b)

file = XDMFFile("poisson.xdmf")
if cpp.common.has_hdf5():
    file.write(w, XDMFFile.Encoding.HDF5)
else:
    file.write(w, XDMFFile.Encoding.ASCII)
