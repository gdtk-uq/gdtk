#!/bin/bash
prep-gas h2-air.inp h2-air-model.lua
gunzip test.dat.gz
onedval test.config test.dat
gzip test.dat
