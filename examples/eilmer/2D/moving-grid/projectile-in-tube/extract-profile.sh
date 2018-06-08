#!/bin/bash

# Let's take the line of data that runs along the bottom of the domain
# in the x-direction.
e4shared --job=projectile-in-tube --post --slice-list="0,:,0,0" --output-file="profile.dat"

# create ideal solution
tclsh piston.tcl

# create graphs
gnuplot plot-profile.gplot

ps2pdf -dEPSCrop pressure-profile.eps
ps2pdf -dEPSCrop temperature-profile.eps
ps2pdf -dEPSCrop velx-profile.eps

ps2pdf -dEPSCrop pressure-time.eps
ps2pdf -dEPSCrop velocity-time.eps
ps2pdf -dEPSCrop position-time.eps
