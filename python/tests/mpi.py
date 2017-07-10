import mpi4py.MPI as pyMPI

print("----")


import dolfin.cpp

print("XXXXX")
#comm = dolfin.cpp.MPI.comm_world


mesh = dolfin.cpp.generation.UnitSquareMesh(2, 2)

tcomm = mesh.mpi_comm()
print(type(tcomm))

#tcomm = dolfin.cpp.MPI.to_comm(tcomm)
#print("****:", type(tcomm))

#dolfin.cpp.MPI.init()

#print("----------")


comm =  pyMPI.COMM_WORLD
#comm =  dolfin.cpp.MPI.comm_world

tcomm = dolfin.cpp.MPI.to_mpi4py_comm(comm)
print(type(tcomm))

#dcomm = test_comm = dolfin.cpp.MPI.to_comm(comm)

#print("----------")
#print("Comm 0:", comm, type(comm))

#print("+++++++++++")
rank = dolfin.cpp.MPI.rank(comm)
print(rank)
rank = dolfin.cpp.MPI.rank(comm)
print(rank)
#print("+++++++++++")


print("^^^^^^^^^^^^")
comm = dolfin.cpp.MPI.comm_world
print("^^^^^^^^^^^^")
print("Comm 1:", comm)

print("$$$$$$$$$$$$$")
size = dolfin.cpp.MPI.size(comm)
print("$$$$$$$$$$$$$")
print("Size:", size)

s = dolfin.cpp.MPI.sum(tcomm, 2.5)
print("Sum 2.5:", s)
