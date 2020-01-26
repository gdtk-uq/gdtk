#!/bin/bash
# prep-chem.sh
cp ${DGD_REPO}/examples/kinetics/hydrogen-ignition-delay/Rogers-Schexnayder-species.inp .
cp ${DGD_REPO}/examples/kinetics/hydrogen-ignition-delay/Rogers-Schexnayder.lua .
sed -i '/^options.H2_O2_only =/s/false/true/' Rogers-Schexnayder-species.inp
sed -i '/^options.H2_O2_only =/s/false/true/' Rogers-Schexnayder.lua
prep-gas Rogers-Schexnayder-species.inp Rogers-Schexnayder-gas-model.lua
prep-chem Rogers-Schexnayder-gas-model.lua Rogers-Schexnayder.lua Rogers-Schexnayder-reac-file.lua
