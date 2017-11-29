#!/bin/bash

prep-gas Evans-Schexnayder-species.inp Evans-Schexnayder-gas-model.lua
prep-chem Evans-Schexnayder-gas-model.lua Evans-Schexnayder-reactions.inp Evans-Schexnayder-chemistry.lua
e4shared --job=lehr --prep

