#!/bin/bash
# prep.sh
DGD_REPO=${DGD_REPO:=${HOME}/gdtk}
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
lrkt-prep --job=ss3
