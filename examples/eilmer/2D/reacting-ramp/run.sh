#! /bin/bash
# run.sh
mpirun -np 8 e4-nk-dist --snapshot-start=last --job=ramp15 # > LOGFILE
#e4-nk-shared --snapshot-start=last --job=ramp15 # > LOGFILE

