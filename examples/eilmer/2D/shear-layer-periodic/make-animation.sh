#! /bin/bash
# make-animation.sh
convert -delay 25 -loop 0 *.png -delay 500 vorticity-field.0100.png shear-layer-vorticity-animation.gif

