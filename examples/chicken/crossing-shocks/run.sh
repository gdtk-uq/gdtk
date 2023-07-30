#! /bin/bash
# run.sh for the crossing-shocks+boundary-layer example
chkn-prep --job=csbli --binary
chkn-run --job=csbli --binary
chkn-post --job=csbli --binary --tindx=all
