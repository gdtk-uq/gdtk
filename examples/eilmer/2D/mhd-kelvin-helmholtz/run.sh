#!/bin/bash
# run.sh

#Runs the Kelvin-Helmholtz Instability with MHD with and without the presence of divergence cleaning
cd KH-MHD-clean

e4shared --prep --job=KH-MHD-clean
if [ "$?" -ne "0" ] ; then
    echo "e4prep-clean ended abnormally."
fi

e4shared --run --job=KH-MHD-clean --verbosity=1 --max-cpus=2
if [ "$?" -ne "0" ] ; then
    echo "e4main-clean ended abnormally."
fi

e4shared --post --job=KH-MHD-clean --tindx-plot=all --vtk-xml
if [ "$?" -ne "0" ] ; then
    echo "e4post-clean ended abnormally."
    exit
fi

cd ..

cd KH-MHD-unclean

e4shared --prep --job=KH-MHD-unclean
if [ "$?" -ne "0" ] ; then
    echo "e4prep-unclean ended abnormally."
fi

e4shared --run --job=KH-MHD-unclean --verbosity=1 --max-cpus=2
if [ "$?" -ne "0" ] ; then
    echo "e4main-unclean ended abnormally."
fi

e4shared --post --job=KH-MHD-unclean --tindx-plot=all --vtk-xml
if [ "$?" -ne "0" ] ; then
    echo "e4post-unclean ended abnormally."
fi

cd ..

e4shared --custom-post --script-file=write_divB_data.lua

gnuplot -e "filename='plotdata.dat'" plot.gnuplot
