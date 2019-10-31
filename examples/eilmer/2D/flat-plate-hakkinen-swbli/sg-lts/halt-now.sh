#! /bin/bash
# halt-now.sh
# Halt the simulation by setting the integer flag
# in the job control file to 1.
# PJ 2019-10-29
cp config/swbli.control config/swbli.control.save-before-halt
sed -e 's/\"halt_now\": 0/\"halt_now\": 1/' \
    config/swbli.control > config/swbli.control.altered
cp config/swbli.control.altered config/swbli.control

