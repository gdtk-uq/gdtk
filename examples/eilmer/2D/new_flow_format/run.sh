#!/bin/bash
# run.sh

# for profiling see https://github.com/prasunanand/gperftools_d:
# $ sudo apt install google-perftools libgoogle-perftools-dev graphviz
# Install go and then install google-pprof.
# $ go get github.com/google/pprof

BINARY=$DGD/bin/e4shared_diagnostic
PPROF=~/go/bin/pprof

prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=ffs


# --- PROFILE RUN ---
CPUPROFILE=./prof.out $BINARY --run --job=ffs --max-cpus=1
$PPROF --pdf $BINARY  ./prof.out > profile.pdf


# -- NORMAL RUN ---
#mpirun -np 3 e4mpi --run --job=ffs
#e4shared --run --job=ffs --verbosity=1 --max-cpus=4
e4shared --post --job=ffs --tindx-plot=all --vtk-xml
e4shared --post --job=ffs --tindx-plot="1,2,3,4,5" --vtk-xml --plotTag="average"
e4shared --post --job=ffs --tindx-plot=last --vtk-xml --plotTag="DFT"
e4shared --post --job=ffs --tindx-plot="1,2,3,4,5" --vtk-xml --plotTag="gradient"

