import mpi4py.MPI as pyMPI

print("----")


import dolfin.cpp

print("XXXXX")
#comm = dolfin.cpp.MPI.comm_world


mesh = dolfin.cpp.generation.UnitSquareMesh(2, 2)

# Uncommenting this messes up later (somethig with memory management)
#tcomm = mesh.mpi_comm()
#print(type(tcomm))

#dolfin.cpp.MPI.init()

#print("----------")


#comm =  pyMPI.COMM_WORLD
comm =  dolfin.cpp.MPI.comm_world
print(comm)

print("A-------------------------------------------------------")
print(type(comm), comm)
tcomm = dolfin.cpp.MPI.to_mpi4py_comm(comm)
print(type(comm), comm)
print("B-------------------------------------------------------")

#print("----------")
print("Comm 0:", comm, type(comm))

print("+++++++++++")
myrank = dolfin.cpp.MPI.size(comm)

print("rank (0a):", comm)
print("rank (0b):", comm, type(tcomm))


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
print(comm)
size = dolfin.cpp.MPI.size(tcomm)
print(comm)
print("$$$$$$$$$$$$$")
print(comm)
print("Size:", size)
print(comm)

print("pre sum", comm)
s = dolfin.cpp.MPI.sum(comm, 2.5)
print("post-sum")
print("Sum 2.5:", s)
print(tcomm)
print(comm)
