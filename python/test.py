import dolfin_test.cpp.geometry
import dolfin_test.cpp.generation
import dolfin_test.cpp.mesh
import dolfin_test.cpp.io
from dolfin_test.function.functionspace import FunctionSpace

# Create some points
p0 = dolfin_test.cpp.geometry.Point(0.0, 0.0, 0.0)
p1 = dolfin_test.cpp.geometry.Point(1.0, 2.0, 4.0)
print(p0)

# Create a mesh
mesh = dolfin_test.cpp.generation.BoxMesh(p0, p1, 2, 6, 9)
print(mesh.num_entities(0))

# Create a mesh
mesh = dolfin_test.cpp.generation.UnitCubeMesh(7, 6, 9)
#print(mesh.num_entities(0))
print(mesh.topology().dim(), mesh.geometry().dim())

# Create a mesh
mesh = dolfin_test.cpp.generation.UnitSquareMesh(6, 9)
print(mesh.num_entities(2))

# Create a mesh
mesh = dolfin_test.cpp.generation.UnitSquareMesh(6, 9, "crossed")
print(mesh.num_entities(2))

# Create a small mesh to check topology and geometry
mesh = dolfin_test.cpp.generation.UnitSquareMesh(2, 2)
print(mesh.topology().dim(), mesh.geometry().dim())
print(mesh.coordinates(), mesh.cells())

Q = FunctionSpace(mesh, "CG", 1)


# Write mesh to file (using two different interfaces)
file = dolfin_test.cpp.io.VTKFile("test.pvd", "ascii")
file.write(mesh)
file << mesh
