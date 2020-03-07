#!/bin/bash
# prep.sh
e4shared --prep --job=bwedge --verbosity=1
e4shared --post --job=bwedge --tindx-plot=0 --vtk-xml
