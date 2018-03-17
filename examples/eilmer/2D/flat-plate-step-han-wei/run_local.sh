#!/bin/bash

e4loadbalance --job=plate --ntasks=8
echo "Start time: ";date>>startlog
mpirun -np 8 e4mpi --job=plate --run
echo "Finish time: ";date>>endlog
