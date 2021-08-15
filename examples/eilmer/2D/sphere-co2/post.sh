#!/bin/bash
# post.sh
#
# 1. flow field pictures
e4shared --post --job=sco2 --tindx-plot=all --vtk-xml --add-vars="mach"
#
# 2. shock location
e4shared --custom-script --script-file=shock-position.lua
