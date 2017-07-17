
import dolfin.cpp as cpp

def _ufl_cell(mesh):
    """Returns the ufl cell of the mesh."""
    import ufl
    gdim = mesh.geometry().dim()
    cellname = mesh.type().description(False)
    return ufl.Cell(cellname, geometric_dimension=gdim)

class UnitIntervalMesh(cpp.generation.UnitIntervalMesh):
    ufl_cell = _ufl_cell

class UnitSquareMesh(cpp.generation.UnitSquareMesh):
    ufl_cell = _ufl_cell

class UnitCubeMesh(cpp.generation.UnitCubeMesh):
    ufl_cell = _ufl_cell


