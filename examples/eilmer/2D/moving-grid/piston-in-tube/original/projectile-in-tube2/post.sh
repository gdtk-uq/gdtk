#!/bin/bash

# We'll ask for all of the flow field snapshots so
# that we can make an animation in Paraview.
# e4shared --job=projectile --post --vtk-xml --tindx-plot=all

# Run the custom post script to check our energy and mass balances
e4shared --custom-script --script-file="balanceCheck.lua"