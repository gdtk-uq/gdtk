#!/bin/bash
e4shared --job=cj-det-wave --post --slice-list="0,:,0,0" --output-file=profile.data
gnuplot plot-comparison.gplot
