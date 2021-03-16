-- prep-flow.lua
-- A place to put helper functions and classes for flow preparation activities.
-- This script is read by the Eilmer4 program at config and flow preparation time,
-- before reading and processing the user's input script.
--
-- Authors: PJ, RJG, Kyle D., Nick G. and Daryl B.
--
print("Loading prep-flow.lua...")

require 'lua_helper'
deepclone = lua_helper.deepclone
checkAllowedNames = lua_helper.checkAllowedNames

require 'blk_conn'
-- Let's pull the symbols out of the blk_conn module
-- and make them global in this namespace
for k,v in pairs(blk_conn) do
   _G[k] = v
end

