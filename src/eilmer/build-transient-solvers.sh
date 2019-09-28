#! /bin/bash
# build-transient-solvers.sh
# Build a debug flavour and a fast flavour of the transient solver.
# 
DGDINST=${HOME}/dgdinst

# We want the debug build to be able to operate with the CEA-gas model
# for the T4 shock tunnel example.
if make DMD=ldc2 FLAVOUR=debug install
then
    echo "Debug build successful"
    if mv ${DGDINST}/bin/e4shared ${DGDINST}/bin/e4shared-debug
    then
        echo "Debug build of shared-memory code installed as e4shared-debug."
    fi
else
    echo "Debug build and install failed."
fi

# The build of the "fast" code will take longer.
make clean
if make DMD=ldc2 FLAVOUR=fast WITH_MPI=1 install
then
    echo "Build and install successful for FLAVOUR=fast."
else
    echo "Build and install failed for FLAVOUR=fast."
fi
