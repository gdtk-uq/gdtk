#! /bin/bash
# build-reduced-transient-solvers.sh
# Build a reduced-capability, fast flavour of the transient solver.
# 
DGDINST=${HOME}/dgdinst

make clean
if make DMD=ldmd2 FLAVOUR=fast WITH_MPI=1 MULTI_SPECIES_GAS=0 MULTI_T_GAS=0 MHD=0 KOMEGA=0 install
then
    echo "Build and install successful for reduced code with FLAVOUR=fast."
else
    echo "Build and install failed for reduced code with FLAVOUR=fast."
fi
