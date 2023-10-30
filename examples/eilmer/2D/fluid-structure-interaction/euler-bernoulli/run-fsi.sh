#!/bin/bash

prep-gas ideal-air.inp ideal-air.lua
e4shared --prep --job=EulerBernoulli
e4shared --run --job=EulerBernoulli
e4shared --post --job=EulerBernoulli --vtk-xml --tindx-plot=all
