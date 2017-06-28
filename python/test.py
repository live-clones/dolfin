#import dolfin_test.generation
import dolfin_test.geometry

#mesh = dolfin_test.generation.UnitSquareMesh(3, 4, "crossed")
#mesh = dolfin_test.generation.UnitCubeMesh(3, 4, 8)

#print(mesh.num_entities(0))
#mesh = dolfin_test.generation.UnitSquareMesh(3, 4, "crossed")

p = dolfin_test.geometry.Point(1.0, 2.0, 4.0)
print(p)
