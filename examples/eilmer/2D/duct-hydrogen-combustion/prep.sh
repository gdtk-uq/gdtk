#!/bin/bash

prep-gas combusting-species.inp h2-o2-n2-9sp.lua
prep-chem h2-o2-n2-9sp.lua Bittker-Scullin.lua h2-o2-n2-9sp-18r.lua
e4shared --job=bittker --prep

