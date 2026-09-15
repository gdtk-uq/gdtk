#!/bin/bash
# run.sh for the Rao thrust-optimised nozzle example.
python3 rao-nozzle.py
gnuplot plot-contour.gnuplot
gnuplot plot-mach.gnuplot
gnuplot plot-theta.gnuplot
