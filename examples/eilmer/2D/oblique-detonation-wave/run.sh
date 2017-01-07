#! /bin/bash
# run.sh
e4shared --prep --job=odw
e4shared --run --job=odw
e4shared --post --job=odw --tindx-plot=all --vtk-xml
