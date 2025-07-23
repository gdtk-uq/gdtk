#!/bin/bash

./ES-test.sh
./RS-test.sh
./Stanford-test.sh
./J92-test.sh

gnuplot plot-ignition-delay.gplot
ps2pdf -dEPSCrop h2-ignition-delay.eps

