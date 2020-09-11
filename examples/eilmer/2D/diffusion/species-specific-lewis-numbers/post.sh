#!/bin/bash

e4shared --job=diffusion --post --slice-list="0,:,0,0" --output-file='profile.data'
gnuplot plot-profile.gplot
