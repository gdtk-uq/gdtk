#!/bin/bash
# make-xt-diagram.sh
l1d4 --xt-data --job=t4-11315 --log10 --var-name=p --tindx-end=9999
echo "Make contour plot of pressure, zooming into shock process."
gnuplot < contour-p-zoom-into-shock.gnuplot
