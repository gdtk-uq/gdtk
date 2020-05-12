#!/bin/bash
# run.sh
prep-gas kiock-gas.inp kiock-gas-model.lua
e4shared --prep --job=kiock
mpirun -np 3 e4mpi --run --job=kiock
