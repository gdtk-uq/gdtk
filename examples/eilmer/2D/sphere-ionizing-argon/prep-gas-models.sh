#!/bin/bash
# prep-gas-models.sh
prep-gas ionizing-argon-model.inp ionizing-argon-model.lua
prep-chem ionizing-argon-model.lua ionizing-argon-reactions.lua \
          ionizing-argon-reactor.lua
