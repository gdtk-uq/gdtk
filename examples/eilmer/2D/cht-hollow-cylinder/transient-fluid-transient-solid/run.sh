# prep gas model
prep-gas ideal-air.inp ideal-air-gas-model.lua

# prep simulation
e4shared --job=cht-cyl --prep

# run MPI version
#mpirun -np 4 e4mpi --job=cht-cyl --run

# run shared memory version
e4shared --job=cht-cyl --run

# post-process simulation
e4shared --job=cht-cyl --post --tindx-plot=all --vtk-xml
