#! /bin/bash
# install-only-debug-transient-solvers.sh
# Build and install only a debug flavour of the transient solver,
# shared-memory variant only.
# 
cd ${DGD_REPO}/src/eilmer
#
# We want the debug build to be able to operate with the CEA-gas model
# for the T4 shock tunnel example.
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
