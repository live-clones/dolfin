import mpi4py.MPI as pyMPI



import dolfin.cpp

#dolfin.cpp.MPI.init()

print("----------")


comm =  pyMPI.COMM_WORLD
#comm = dolfin.cpp.MPI.comm_world

dcomm = test_comm = dolfin.cpp.MPI.to_comm(comm)

#print("----------")
#print("Comm 0:", comm)

print("+++++++++++")
rank = dolfin.cpp.MPI.rank(comm)
print(rank)
rank = dolfin.cpp.MPI.rank(dcomm)
print(rank)
print("+++++++++++")


"""
print("^^^^^^^^^^^^")
comm = dolfin.cpp.MPI.comm_world
print("^^^^^^^^^^^^")
print("Comm 1:", comm)

print("$$$$$$$$$$$$$")
size = dolfin.cpp.MPI.size(comm)
print("$$$$$$$$$$$$$")
print("Size:", size)

s = dolfin.cpp.MPI.sum(comm, 2.5)
print("Sum 2.5:", s)
"""
