#!/bin/bash
# We'll ask for all of the flow field snapshots so
# that we can make an animation in Paraview.
e4shared --job=piston --post --vtk-xml --tindx-plot=all
