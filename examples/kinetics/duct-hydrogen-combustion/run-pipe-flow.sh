#!/bin/bash
# run-pipe-flow.sh
prep-gas combusting-species.inp h2-o2-n2-9sp.lua
prep-chem h2-o2-n2-9sp.lua Bittker-Scullin.lua h2-o2-n2-9sp-18r.lua
gas-calc reacting_pipe_flow.lua > reacting-pipe-flow.transcript
# Make any non-data lines into comments.
sed -i '/^[a-zA-Z]/s/^/# /' reacting-pipe-flow.transcript
gnuplot plot-profiles.gplot
