import dolfin.cpp
import numpy as np


def create_meshview(mesh_function, value):
    """Create a MeshView of the entities of a MeshFunction marked with a given value

    :param mesh_function: The mesh function
    :param value: The marker value
    :return: The MeshView
    """
    current_log_level = dolfin.cpp.log.get_log_level()
    mesh = mesh_function.mesh()
    dim = mesh_function.dim()
    # Codim 0: Normal construction
    if dim == mesh.topology().dim():
        mv = dolfin.cpp.mesh.MeshView.create(mesh_function, value)

        # Create bounding-box tree on all processes to avoid hanging at assembly
        dolfin.cpp.log.set_log_level(dolfin.cpp.log.LogLevel.WARNING)
        mv.bounding_box_tree()
        dolfin.cpp.log.set_log_level(current_log_level)
        return mv
    else:
        assert dim == mesh.topology().dim() - 1, "MeshViews of codim > 1 not supported"

    # Codim 1: Facets does not have ownership, and is assigned to the lowest
    # rank that has them
    shared_facets = mesh.topology().shared_entities(dim)
    self_rank = mesh.mpi_comm().rank
    mesh_view_mf = dolfin.cpp.mesh.MeshFunctionSizet(mesh, dim, 0)
    for i in range(len(mesh_function.array())):
        shared = shared_facets.get(i, None)
        if shared is not None:
            shared_procs = np.asarray(list(shared), dtype=np.int64)
            if min(shared_procs) > self_rank:
                mesh_view_mf.array()[i] = mesh_function.array()[i]
        else:
            mesh_view_mf.array()[i] = mesh_function.array()[i]

    mv = dolfin.cpp.mesh.MeshView.create(mesh_view_mf, 1)

    # Create bounding-box tree on all processes to avoid hanging at assembly
    dolfin.cpp.log.set_log_level(dolfin.cpp.log.LogLevel.WARNING)
    mv.bounding_box_tree()
    dolfin.cpp.log.set_log_level(current_log_level)

    return mv
