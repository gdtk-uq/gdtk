#!/bin/bash

e4shared --job=sharma --post --vtk-xml --tindx-plot=all

e4shared --job=sharma --post --tindx-plot=last \
         --slice-list="0,:,0,0;1,:,0,0" --output-file=axis-Blackman.data
