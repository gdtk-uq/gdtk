#!/bin/bash
# make-xt-diagram-matplotlib.sh
l1d4 --xt-data-json --job=he-air --tindx-end=9999 --millisec
echo "Make contour plot of pressure over full facility."
l1d4-xt-viewer --variable=p --take-log --v-range=3.0:0.1:6.6 \
               --x-range=-3.0:1.0:8.5 *.json
