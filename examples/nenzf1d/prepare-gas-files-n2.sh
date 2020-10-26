#! /bin/bash
# prepare-gas-files-n2.sh
cp ${DGD_REPO}/examples/kinetics/fixed-volume-reactor-n2/cea-n2-gas-model.lua .
cp ${DGD_REPO}/examples/kinetics/fixed-volume-reactor-n2/nitrogen-2sp.inp .
cp ${DGD_REPO}/examples/kinetics/fixed-volume-reactor-n2/nitrogen-2sp-2r.lua .
prep-gas nitrogen-2sp.inp nitrogen-2sp-1T.lua
prep-chem nitrogen-2sp-1T.lua nitrogen-2sp-2r.lua nitrogen-2sp-1T-reactions.lua
