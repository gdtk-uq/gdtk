#!/usr/bin/env python
from __future__ import print_function

import os

cmd = "prep-gas ideal-air.inp ideal-air.lua"
os.system(cmd)

cmd = "e4shared --job=reference --prep"
os.system(cmd)

cmd = "e4shared --job=reference --run"
os.system(cmd)

probe_str = ""
for j in range(4):
    for k in range(4):
        probe_str += "%f,%f,%f;" % (1.125, 0.125+0.25*j, 0.125+0.25*k)
cmd = 'e4shared --job=reference --post --tindx-plot=last --probe="%s" --output-file=ref.out' % probe_str[:-1]
os.system(cmd)

print("Generated plane of reference data: ref.out")



