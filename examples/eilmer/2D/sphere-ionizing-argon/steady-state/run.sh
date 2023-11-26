#!/bin/bash
# run.sh
#e4-nk-shared --max-cpus=1 --job=sphere --snapshot-start=last
mpirun -np 12 e4-nk-dist --job=sphere --snapshot-start=last
