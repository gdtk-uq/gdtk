#!/bin/bash

# Extract the solution data over whole flow domain and reformat.
e4shared --post --job=plate --tindx-plot=last --vtk-xml --add-vars="mach,total-h"

