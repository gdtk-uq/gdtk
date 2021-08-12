#!/bin/bash
# run.sh

# for profiling see https://github.com/prasunanand/gperftools_d:
# $ sudo apt install google-perftools libgoogle-perftools-dev graphviz
# Install go and then install google-pprof.
# $ go get github.com/google/pprof

BINARY=$DGD/bin/e4shared_diagnostic
PPROF=~/go/bin/pprof

FORMAT=rawbinary

if  true ; then
    cp ffs.lua ffs_run.lua
    
    if [[ ${FORMAT} == *"eilmer4"* ]]; then
        printf "\nconfig.new_flow_format = true" >> ffs_run.lua
    fi
    
    printf "\nconfig.flow_format = \"%s\"" ${FORMAT} >>  ffs_run.lua

    prep-gas ideal-air.inp ideal-air-gas-model.lua
    e4shared --prep --job=ffs_run


    # --- PROFILE RUN ---
    CPUPROFILE=./prof.out $BINARY --run --job=ffs_run --max-cpus=1

fi

$PPROF --pdf --nodecount=100  $BINARY  ./prof.out > profile_${FORMAT}.pdf
