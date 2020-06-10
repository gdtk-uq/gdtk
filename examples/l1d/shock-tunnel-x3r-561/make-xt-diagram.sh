#!/bin/bash
# make-xt-diagram.sh
l1d4 --xt-data --job=x3r-561 --log10 --var-name=p --tindx-end=9999
echo "Make contour plot of pressure over full facility."
gnuplot < contour-p.gnuplot

l1d4 --xt-data --job=x3r-561 --var-name=T --tindx-end=9999
echo "Make contour plot of temperature, over full facility."
gnuplot < contour-T.gnuplot
