#! /bin/bash
# build-complex-solvers.sh
#
# Trial build of Kyle's complex-numbers flavour of the solvers.

if make DMD=ldmd2  WITH_COMPLEX_NUMBERS=1 WITH_SSC=1 WITH_SSS=1
then
    echo "Build successful for complex flavour."
else
    echo "Build failed."
fi
