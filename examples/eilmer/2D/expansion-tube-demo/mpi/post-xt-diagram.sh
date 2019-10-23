#!/bin/bash
# post-xt-diagram.sh
xtdata.rb --job=exptube --log10 --xcolumn=1 --vcolumn=9 --output=xtlogp.data
gnuplot < contour.gnuplot
