#!/bin/bash
# run.sh
# e4shared --run --job=ss3 --verbosity=1 --max-cpus=4
mpirun -np 4 e4mpi --run --job=ss3 --verbosity=1
e4shared --post --job=ss3 --tindx-plot=all --vtk-xml \
	 --add-vars="mach,pitot,total-p,total-h"
