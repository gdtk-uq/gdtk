#!/bin/bash

gunzip test.dat.gz
onedval test.config test.dat 
gzip test.dat

