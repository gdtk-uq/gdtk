#!/bin/bash
cp ${DGD_REPO}/examples/kinetics/nitrogen-dissociation-2T/nitrogen-2sp-2T.inp .
cp ${DGD_REPO}/examples/kinetics/nitrogen-dissociation-2T/Park93-nitrogen-dissociation-2T.lua .
cp ${DGD_REPO}/examples/kinetics/nitrogen-dissociation-2T/N2-energy-exchange.lua .

prep-gas nitrogen-2sp-2T.inp nitrogen-gas-model.lua
prep-chem nitrogen-gas-model.lua Park93-nitrogen-dissociation-2T.lua N2-diss-2T.chem
e4shared --job=n90 --prep
