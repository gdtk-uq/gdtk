#!/bin/bash

# prep
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cavity
