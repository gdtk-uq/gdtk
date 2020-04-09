#! /bin/bash
# install-loadable-module.sh
# Build and install a dynamically-loadable module of gas and geometry functions.
# There are Python and Ruby packages that use/import this module.
# 
cd ${DGD_REPO}/src/gas
#
make clean
if make DMD=dmd install
then
    echo "Build and install successful for loadable module."
else
    echo "Build and install failed for loadable module."
fi
