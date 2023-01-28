#!/bin/bash
# display-xt-diagram.sh
l1d4 --xt-data-json --job=he-air-1T --tindx-end=9999
l1d4-xt-viewer --variable p --take-log --v-range 2.4:0.1:6.0 *.json
