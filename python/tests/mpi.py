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


#comm =  pyMPI.COMM_WORLD
comm =  dolfin.cpp.MPI.comm_world

print("A-------------------------------------------------------")
#print(type(comm))
tcomm = dolfin.cpp.MPI.to_mpi4py_comm(comm)
print(type(tcomm))
print("B-------------------------------------------------------")

#print("----------")
print("Comm 0:", comm, type(tcomm))

#print("+++++++++++")
myrank = dolfin.cpp.MPI.size(tcomm)

print("rank (0a):", myrank)
print("rank (0b):", type(tcomm), tcomm)


#print("C-------------------------------------------------------")
#rank = dolfin.cpp.MPI.rank(tcomm)
#print("rank (1a):", rank)
#print("rank (1b):", type(tcomm))

#print("rank (1):", rank, type(tcomm))
#print("+++++++++++")


#print("^^^^^^^^^^^^")
#comm = dolfin.cpp.MPI.comm_world
#print("^^^^^^^^^^^^")
#print("Comm 1:", comm)

print("$$$$$$$$$$$$$")
size = dolfin.cpp.MPI.size(comm)
print("$$$$$$$$$$$$$")
print("Size:", size)

print(type(tcomm))
s = dolfin.cpp.MPI.sum(tcomm, 2.5)
print("Sum 2.5:", s)
print(tcomm)
