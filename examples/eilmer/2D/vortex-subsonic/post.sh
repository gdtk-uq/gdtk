#!/bin/bash
# post.sh
e4shared --post --job=vortex --vtk-xml --tindx-plot=all \
         --add-vars="mach,total-h,entropy"
