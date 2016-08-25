#!/bin/bash
e4shared --custom-post --script-file=reacting_pipe_flow.lua > reacting-pipe-flow.transcript
sed -i '/^[a-zA-Z]/s/^/# /' reacting-pipe-flow.transcript
