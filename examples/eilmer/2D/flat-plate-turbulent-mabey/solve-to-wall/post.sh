#!/bin/bash

e4shared --post --job=mabey --output-file="mabey-x-368mm.dat" \
  --slice-list="2,28,:,0;3,28,:,0" --add-vars="mach,pitot,total-p,total-h"

e4shared --post --job=mabey --output-file="mabey-y-wall.dat" \
  --slice-list="1,:,42,0;3,:,42,0" --add-vars="mach,pitot,total-p,total-h"
