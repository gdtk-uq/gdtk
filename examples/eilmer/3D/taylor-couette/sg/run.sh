#!/bin/bash
sleep 10 && gnuplot live-residuals.gplot &
e4sss --job=low --max-cpus=4

