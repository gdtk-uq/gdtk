#! /bin/bash
# run.sh for the short flat-plate boundary-layer example
chkn-prep --job=short --binary
chkn-run --job=short --binary
chkn-post --job=short --binary --tindx=all
chkn-post --job=short --binary --probe=1.0,0.0,0.0
chkn-post --job=short --binary --slice="0,0,0,84,1,:"
