#!/bin/bash

# We'll ask for all of the flow field snapshots so
# that we can make an animation in Paraview.

e4shared --job=cube --post --vtk-xml --tindx-plot=all --add-vars="mach,pitot,total-p,total-h" --output-file="cube_surfaces" --surface-list="0,5;3,5;6,5;9,5;12,5;15,5;18,5;21,5;24,5;81,5;84,5;87,5;90,5;93,5;96,5;99,5;102,5;105,5;108,5;111,5;114,5;117,5;120,5;123,5;126,5;129,5;132,5;29,4;32,4;35,4;38,4;41,4;44,4;47,4;50,4;53,4;56,4;59,4;62,4;65,4;68,4;71,4;74,4;77,4;80,4;137,4;140,4;143,4;146,4;149,4;152,4;155,4;158,4;161,4"
