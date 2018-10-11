#!/bin/bash

prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --job=stokes-2nd-flow --prep
