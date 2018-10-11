#!/bin/bash

# Let's take the line of data that runs along the bottom of the domain
# in the x-direction.
e4shared --job=stokes-2nd-flow --post --slice-list="0,0,:,0" --tindx-plot=0  --output-file="profile0.dat"
e4shared --job=stokes-2nd-flow --post --slice-list="0,0,:,0" --tindx-plot=10 --output-file="profile1.dat"
e4shared --job=stokes-2nd-flow --post --slice-list="0,0,:,0" --tindx-plot=20 --output-file="profile2.dat"
e4shared --job=stokes-2nd-flow --post --slice-list="0,0,:,0" --tindx-plot=30 --output-file="profile3.dat"
e4shared --job=stokes-2nd-flow --post --slice-list="0,0,:,0" --tindx-plot=40 --output-file="profile4.dat"
e4shared --job=stokes-2nd-flow --post --slice-list="0,0,:,0" --tindx-plot=50 --output-file="profile5.dat"

# create ideal solution
tclsh stokesflow.tcl

# create graphs
gnuplot plot-profile.gplot

ps2pdf -dEPSCrop velocity-profiles.eps

