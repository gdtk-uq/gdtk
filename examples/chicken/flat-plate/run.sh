#! /bin/bash
# run.sh for the flat-plate boundary-layer example
chkn-prep --job=plate --binary
chkn-run --job=plate --binary
chkn-post --job=plate --binary --tindx=all
chkn-post --job=plate --binary --probe=1.0,0.0,0.0
chkn-post --job=plate --binary --slice="0,0,0,209,1,:"
