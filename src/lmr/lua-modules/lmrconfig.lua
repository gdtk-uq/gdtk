-- lmr_config.lua
--
-- A module for helper functions to query and interact with
-- Eilmer program config. Here, we are referring to static
-- configuration of the software itself, not the configuration
-- related to individual simulations.
--
-- Authors: RJG, PJ, KAD, NNG
--

local lmr_config = {}

local json = require 'json'

function lmrConfigAsTable()
   local lmrCfgFile = os.getenv("DGD") .. "/etc/lmr.cfg"
   local f = assert(io.open(lmrCfgFile, "r"))
   local jsonStr = f:read("*a")
   f:close()
   local jsonData = json.parse(jsonStr)
   return jsonData
end

local lmrCfg = lmrConfigAsTable()
lmr_config.lmrCfg = lmrCfg

function lmr_config.simulationConfigFilename()
   return lmrCfg["config-directory"] .. "/" .. lmrCfg["config-filename"]
end

function lmr_config.simulationControlFilename()
   return lmrCfg["config-directory"] .. "/" .. lmrCfg["control-filename"]
end

function lmr_config.blockListFilename()
   return lmrCfg["config-directory"] .. "/" .. lmrCfg["block-list-filename"]
end

function lmr_config.nkConfigFilename()
   return lmrCfg["config-directory"] .. "/" .. lmrCfg["newton-krylov-config-filename"]
end

function lmr_config.gridMetadataFilename(id)
   -- When id is supplied give individual block
   if id then
      return lmrCfg["grid-directory"] .. "/" .. string.format(lmrCfg["block-filename-format"], id) .. "." .. lmrCfg["grid-metadata-filename"]
   end
   -- else, return global metadataname   
   return lmrCfg["grid-directory"] .. "/" .. lmrCfg["grid-metadata-filename"]
end

function lmr_config.gridFilenameWithoutExt(id)
   return lmrCfg["grid-directory"] .. "/" .. string.format(lmrCfg["block-filename-format"], id)
end

function lmr_config.steadyFlowDirectory(snapshot)
   dname = lmrCfg["snapshot-directory"]
   dname = dname .. "/"
   dname = dname .. string.format(lmrCfg["snapshot-index-format"], snapshot)
   dname = dname .. "/"
   dname = dname .. lmrCfg["flow-directory"]
   return dname
end

function lmr_config.steadyFlowFilename(snapshot, blkId)
   fname = lmr_config.steadyFlowDirectory(snapshot)
   fname = fname .. "/"
   fname = fname .. string.format(lmrCfg["block-filename-format"], blkId)
   fname = fname .. "."
   fname = fname .. lmrCfg["zip-extension"]
   return fname
end

return lmr_config

