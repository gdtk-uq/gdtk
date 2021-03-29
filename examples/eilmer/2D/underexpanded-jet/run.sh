#!/bin/bash
# run.sh
e4shared --run --job=uej --verbosity=1 --max-cpus=8
-- mpirun -oversubscribe -np 7 e4mpi --run --job=uej --verbosity=1
