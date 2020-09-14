#!/bin/bash
# post-profile.sh
e4shared --post --job=inlet --tindx-plot=last \
         --add-vars="mach,total-p" \
         --slice-list="1,$,:,0" \
         --output-file="profile.data"
