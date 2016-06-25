#!/usr/bin/python
from dolfin import *

PETScOptions.set("ksp_view");
PETScOptions.set("ksp_monitor_true_residual");
PETScOptions.set("pc_type", "fieldsplit");
PETScOptions.set("pc_fieldsplit_type", "additive");
PETScOptions.set("fieldsplit_0_pc_factor_mat_solver_package", "superlu_dist");
PETScOptions.set("fieldsplit_0_ksp_type", "preonly");
PETScOptions.set("fieldsplit_1_pc_factor_mat_solver_package", "superlu_dist");
PETScOptions.set("fieldsplit_1_ksp_type", "preonly");
PETScOptions.set("fieldsplit_0_pc_type", "lu");
PETScOptions.set("fieldsplit_1_pc_type", "lu");

# Sub domain for top and bottom
def TopBottom(x, on_boundary):
    return abs(1.0 - x[1]) < DOLFIN_EPS or abs(x[1]) < DOLFIN_EPS

def LeftEdge(x, on_boundary):
    return x[0] < DOLFIN_EPS

mesh = UnitSquareMesh(50, 50);

# Create function spaces
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, P2)
Q = FunctionSpace(mesh, P1)

u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Define forms for each block
a00 = inner(grad(u), grad(v))*dx
a01 = - div(v)*p*dx
a10 = - div(u)*q*dx
# Only need a11 if we have BCs in Q space
a11 = Constant(0.0)*p*q*dx
p11 = p*q*dx

f = Constant((0.0, 1.0))
g = Constant(0.0)
L0 = dot(f, v)*dx
L1 = dot(g, q)*dx

# Velocity BC
flow_velocity = Constant((0.0, 0.0));
bc0 = DirichletBC(V, flow_velocity, TopBottom)

# Pressure BC
zero = Constant(12.0)
bc1 = DirichletBC(Q, zero, LeftEdge)

# Assemble all blocks with BCs
bcsV = [bc0]
bcsQ = []

assemblerA = SystemAssembler([a00, a01, a10, None], [L0, L1], [bcsV, bcsQ])

A00 = PETScMatrix()
A01 = PETScMatrix()
A10 = PETScMatrix()

b0 = PETScVector()
b1 = PETScVector()

assemblerA.assemble([A00, A01, A10, None], [b0, b1])

P11 = PETScMatrix()
assemblerP = SystemAssembler(p11, L1, bcsQ)
assemblerP.assemble(P11)

print "A00 = ", A00.size(0), "x" , A00.size(1), A00.norm("frobenius")
print "A01 = ", A01.size(0), "x" , A01.size(1), A01.norm("frobenius")
print "A10 = ", A10.size(0), "x" , A10.size(1), A10.norm("frobenius")
print "P11 = ", P11.size(0), "x" , P11.size(1), P11.norm("frobenius")

# Combine matrices
AA = PETScNestMatrix([A00, A01, A10, None])
PP = PETScNestMatrix([A00, None, None, P11])

# Create solution vectors
u = Function(V)
p = Function(Q)
# Map RHS and solution into block space vectors
x = Vector()
b = Vector()
AA.init_vectors(x, [u.vector(), p.vector()])
AA.init_vectors(b, [b0, b1])

solver = PETScKrylovSolver("minres")
solver.set_from_options()
solver.set_operators(AA, PP)


PETScPreconditioner.set_fieldsplit(solver, AA, ["0", "1"]);

solver.solve(x, b)

print "u.norm = ", u.vector().norm("l2")
