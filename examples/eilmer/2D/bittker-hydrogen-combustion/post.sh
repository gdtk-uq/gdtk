#!/bin/bash

e4shared --job=bittker --post --output-file="profile.data" --slice-list=":,:,0,0" --tindx-plot=last
gnuplot plot-profiles.gplot

