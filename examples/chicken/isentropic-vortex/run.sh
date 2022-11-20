#! /bin/bash
# run.sh for the isentropic vortex example
chkn-prep --job=vortex --binary
chkn-run --job=vortex --binary
chkn-post --job=vortex --binary --tindx=all
