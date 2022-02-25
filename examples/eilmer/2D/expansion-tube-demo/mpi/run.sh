#!/bin/bash
# run.sh
mpirun -np 8 --oversubscribe e4mpi --run --job=exptube --verbosity=1
