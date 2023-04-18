#! /bin/bash
# install-solvers.sh
# Build and install the transient and steady-state solvers
# using the same set of build options as Nick's CI bot.
#   e4shared
#   e4monitor
#   e4zshared
#   e4-nk-shared
#   e4-nk-shared-real
#   e4mpi
#   e4loadbalance
#   e4zmpi
#   e4-nk-dist
#   e4-nk-dist-real
#   e4shared-debug
#
cd ${DGD_REPO}/src/eilmer
#
if make DMD=ldc2 \
        WITH_MPI=1 \
        WITH_COMPLEX_NUMBERS=1 \
        WITH_NK=1 \
        WITH_E4DEBUG=1 \
        FLAVOUR=fast \
        install
then
    echo "Build and install successful for FLAVOUR=fast."
else
    echo "Build and install failed for FLAVOUR=fast."
fi
