#!/bin/bash
# prep.sh
DGD_REPO=${DGD_REPO:=${HOME}/dgd}
cp ${DGD_REPO}/src/gas/sample-data/cea-lut-air-version-test.lua ./cea-lut-air.lua
e4shared --prep --job=ss3
