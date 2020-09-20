#! /bin/bash
# post.sh
e4shared --post --job=tbo --vtk-xml \
         --tindx-plot=all --add-vars="total-p,mach"
