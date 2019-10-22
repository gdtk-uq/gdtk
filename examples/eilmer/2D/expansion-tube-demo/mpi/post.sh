#!/bin/bash
# post.sh
e4shared --post --job=exptube --vtk-xml --tindx-plot=all \
	 --add-vars="mach,pitot,total-p,total-h"
