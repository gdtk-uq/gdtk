#!/bin/bash

./ES-test.sh
./RS-test.sh
./Stanford-test.sh

gnuplot plot-ignition-delay.gplot
ps2pdf -dEPSCrop h2-ignition-delay.eps

