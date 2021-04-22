#!/bin/bash
# make-xt-diagram.sh
l1d4 --xt-data --job=x2_condition1 --log10 --var-name=p --tindx-end=9999
gnuplot < contour-p.gnuplot
