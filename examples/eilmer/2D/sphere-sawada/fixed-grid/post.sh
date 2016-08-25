#! /bin/sh
# post.sh
e4shared --post --job=ss3 --tindx-plot=last --slice-list="0,:,0,:" \
         --output-file=ss3_stag_line.data
awk -f locate_shock.awk ss3_stag_line.data > ss3.result
