#!/bin/bash
# run.sh
prep-gas idealN2.inp ideal-N2-gas-model.lua
e4shared --prep --job=cyl-flare
e4shared --run --job=cyl-flare --verbosity=1
e4shared --post --job=cyl-flare --tindx-plot=last --extract-line="0.000112445,0.0325093,0,0.101325,0.0325369,0,1000;0.101987,0.0327078,0,0.219704,0.100676,0,1000" --output-file="wallData" --vtk-xml --add-vars="mach,pitot,total-p,total-h"

