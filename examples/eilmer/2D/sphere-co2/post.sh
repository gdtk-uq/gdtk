#!/bin/bash
# post.sh
#
# 1. flow field pictures
e4shared --post --job=sco2 --tindx-plot=all --vtk-xml --add-vars="mach"
#
# 2. shock location
e4shared --custom-script --script-file=shock-position.lua
#
# 3. stagnation-line properties (first line of cell centres)
e4shared --post --job=sco2 --tindx-plot=last --output-file="stag-profile.data" \
         --slice-list="0,:,0,0;2,:,0,0" --add-vars="mach,pitot"
gnuplot plot-stag-line-properties.gplot
