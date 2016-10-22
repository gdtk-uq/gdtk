#!/bin/bash
# post.sh

# Create a VTK plot file of the steady full flow field.
e4shared --post --job=cyl --tindx-plot=last --vtk-xml

# Pull out the cylinder surfaces.
e4shared --post --job=cyl --tindx-plot=last --output-file=cylinder \
    --add-vars="mach" --surface-list="0,east;1,east;3,bottom"

# Now pull out some block surfaces that show cross-sections of the flow field.
e4shared --post --job=cyl --tindx-plot=last --output-file=interior \
    --add-vars="mach" \
    --surface-list="0,bottom;1,bottom;0,north;1,north;2,north;3,north;0,south;1,south;2,south;3,south;3,east"

# Stagnation-line flow data
e4shared --post --job=cyl --tindx-plot=last --output-file=stagnation-line.data \
    --add-vars="mach" --slice-list="0,:,0,0" \
