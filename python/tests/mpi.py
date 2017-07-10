import dolfin.cpp

dolfin.cpp.MPI.init()

print("----------")
comm = dolfin.cpp.MPI.comm_world
print("----------")
print("Comm 0:", comm)

print("+++++++++++")
rank = dolfin.cpp.MPI.rank(comm)
print("+++++++++++")
print("Comm 1:", type(comm), comm)

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
