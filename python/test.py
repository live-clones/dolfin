import dolfin_test.generation
import dolfin_test.geometry
import dolfin_test.mesh
import dolfin_test.io

# Create some points
p0 = dolfin_test.geometry.Point(0.0, 0.0, 0.0)
p1 = dolfin_test.geometry.Point(1.0, 2.0, 4.0)
print(p0)

# Create a mesh
mesh = dolfin_test.generation.BoxMesh(p0, p1, 2, 6, 9)
print(mesh.num_entities(0))

# Create a mesh
mesh = dolfin_test.generation.UnitCubeMesh(7, 6, 9)
print(mesh.num_entities(0))

# Create a mesh
mesh = dolfin_test.generation.UnitSquareMesh(6, 9)
print(mesh.num_entities(2))

# Create a mesh
mesh = dolfin_test.generation.UnitSquareMesh(6, 9, "crossed")
print(mesh.num_entities(2))

# Write mesh to file
file = dolfin_test.io.VTKFile("test.pvd", "ascii")
file.write(mesh)
file << mesh
