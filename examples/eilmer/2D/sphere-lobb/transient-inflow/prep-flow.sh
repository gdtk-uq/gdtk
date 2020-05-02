#!/bin/bash
# prep-flow.sh
e4shared --prep --job=lobb --verbosity=1
e4shared --post --job=lobb --vtk-xml --tindx-plot=0
