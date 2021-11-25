#!/bin/bash

# We'll ask for all of the flow field snapshots so
# that we can make an animation in Paraview.
e4shared --job=cube --post --vtk-xml --tindx-plot=all --add-vars="mach,pitot,total-p,total-h"
