#!/bin/bash
# prep-chem.sh
prep-gas species.inp gas-model.lua
prep-chem gas-model.lua reaction_mechanism.lua reac-file.lua
e4shared --job=ramp15 --prep
e4shared --job=ramp15 --post --tindx-plot=last --vtk-xml --add-vars="mach"
