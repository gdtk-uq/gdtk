#! /bin/bash
# install-transient-solvers.sh
# Build and install a debug flavour of the shared-memory solver and
# build and install a fast flavour of the transient solver,
# shared-memory and MPI variants.
#
cd ${DGD_REPO}/src/eilmer
#
# We want the debug build to be able to operate with the CEA-gas model
# for the T4 shock tunnel example.
# This build should be relatively quick.
if make DMD=ldc2 FLAVOUR=debug install
then
    echo "Debug build successful"
    if mv ${DGD}/bin/e4shared ${DGD}/bin/e4shared-debug
    then
        echo "Debug build of shared-memory code installed as e4shared-debug."
    fi
else
    echo "Debug build and install failed."
fi
#
# The build of the "fast" code will take longer.
make clean
if make DMD=ldc2 FLAVOUR=fast WITH_MPI=1 install
then
    echo "Build and install successful for FLAVOUR=fast."
else
    echo "Build and install failed for FLAVOUR=fast."
fi
