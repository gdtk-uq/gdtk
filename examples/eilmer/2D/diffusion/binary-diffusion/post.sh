#!/bin/bash

e4shared --job=bd --post --slice-list="0,:,0,0" --output-file='profile.data'
gnuplot plot-profile.gplot
