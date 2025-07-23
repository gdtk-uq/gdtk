#!/bin/bash
cp $DGD_REPO/examples/kinetics/hydrogen-ignition-delay/Jachimowski-1992-species.inp ./
cp $DGD_REPO/examples/kinetics/hydrogen-ignition-delay/Jachimowski-1992.lua ./

prep-gas Jachimowski-1992-species.inp gm.lua
prep-chem gm.lua Jachimowski-1992.lua rr.lua

slf hydrogen.yaml

python3 ../scripts/residuals.py hydrogen.log
python3 ../scripts/solution.py hydrogen.sol
