#! /bin/bash
# plot.sh
e4shared --post --job=swbli --tindx-plot=last --add-vars="mach" --output-file=bl.data \
    --extract-line="5.5e-05,1.8e-05,0,0.089,1.8e-05,0,1000"
awk -f compute-shear.awk bl.data > shear.data
gnuplot pressure.gnuplot
gnuplot shear.gnuplot
echo "At this point, we should have pictures to view"


