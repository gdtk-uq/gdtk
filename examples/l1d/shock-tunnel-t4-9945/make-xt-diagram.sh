#!/bin/bash
# make-xt-diagram.sh
l1d4 --xt-data --job=t4-9945 --log10 --var-name=p --tindx-end=9999
echo "Make contour plot of pressure over full facility."
gnuplot < contour-p.gnuplot
echo "Make contour plot of pressure, zooming into shock process."
gnuplot < contour-p-zoom-into-shock.gnuplot

l1d4 --xt-data --job=t4-9945 --var-name=T --tindx-end=9999
echo "Make contour plot of temperature, zooming into shock process."
gnuplot < contour-T-zoom-into-shock.gnuplot
