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