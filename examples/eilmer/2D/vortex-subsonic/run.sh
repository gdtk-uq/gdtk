#!/bin/bash
# run.sh
# e4shared --run --job=vortex --verbosity=1

# Test OpenMPI version. If 4.0 or later, we need to explicitly ask for hwthreads.
OMPI_VER=$(ompi_info --parsable | grep ompi:version:full | awk -F: '{print $4}' | awk -F. '{print $1}')

if [ $OMPI_VER -ge 4 ]
then
    mpirun --use-hwthread-cpus -np 4 e4mpi --run --job=vortex --verbosity=1
else
    mpirun -np 4 e4mpi --run --job=vortex --verbosity=1
fi
