#!/bin/bash
# prep-chem.sh
cp ${DGD_REPO}/examples/kinetics/hydrogen-ignition-delay/Rogers-Schexnayder-species.inp .
cp ${DGD_REPO}/examples/kinetics/hydrogen-ignition-delay/Rogers-Schexnayder.lua .
prep-gas Rogers-Schexnayder-species.inp Rogers-Schexnayder-gas-model.lua
prep-chem Rogers-Schexnayder-gas-model.lua Rogers-Schexnayder.lua Rogers-Schexnayder-reac-file.lua
