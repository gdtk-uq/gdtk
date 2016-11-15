#!/bin/bash
sleep 1 && gnuplot live-residuals.gplot &
e4sss --job=wedge --max-CPUs=2

