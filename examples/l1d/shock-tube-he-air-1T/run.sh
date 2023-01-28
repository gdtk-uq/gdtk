# run.sh
# Tamara's 1-D shock tube exercise, helium driving 1T air.
#
l1d4-prep --job=he-air-1T
l1d4 --run-simulation --job=he-air-1T
