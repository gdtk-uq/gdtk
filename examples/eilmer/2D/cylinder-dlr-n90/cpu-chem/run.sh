#!/bin/bash
# run.sh
e4shared --prep --job=n90
mpirun -np 4 e4mpi --run --job=n90 --verbosity=1
e4shared --post --job=n90 --tindx-plot=all --vtk-xml
