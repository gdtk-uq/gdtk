#!/bin/bash
# make-xt-diagram.sh
l1d4 --xt-data --job=dn2 --log10 --var-name=p --tindx-end=9999
echo "Make contour plot of pressure over full facility."
gnuplot < contour-p.gnuplot
