#! /bin/bash
# install-reduced-transient-solvers-debug.sh
# Build and install a reduced-capability, debug flavour of the transient solver,
# both shared-memory and MPI variants.
# 
cd ${DGD_REPO}/src/eilmer
#
make clean
if make DMD=ldmd2 FLAVOUR=debug WITH_MPI=1 MULTI_SPECIES_GAS=0 MULTI_T_GAS=0 MHD=0 TURBULENCE=0 install
then
    echo "Build and install successful for reduced code with FLAVOUR=debug."
else
    echo "Build and install failed for reduced code with FLAVOUR=debug."
fi
