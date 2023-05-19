#! /bin/bash
# run.sh for the injector example
chkn-prep --job=inj --binary
chkn-run --job=inj --binary
chkn-post --job=inj --binary --tindx=all
