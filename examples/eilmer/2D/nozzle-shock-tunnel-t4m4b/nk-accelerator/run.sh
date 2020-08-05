#!/bin/bash
# run.sh
mpirun -np 4 e4-nk-dist --snapshot-start=last --job=t4m4b --verbosity=1
#e4-nk-shared --snapshot-start=last --job=t4m4b --verbosity=1

