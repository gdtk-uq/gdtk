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

local globalconfig = require 'globalconfig'
config = globalconfig.config

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

function lmr_config.mpimapFilename()
   return lmrCfg["config-directory"] .. "/" .. lmrCfg["mpimap-filename"]
end


function lmr_config.nkConfigFilename()
   return lmrCfg["config-directory"] .. "/" .. lmrCfg["newton-krylov-config-filename"]
end

function lmr_config.gridDirectory()
   return lmrCfg["simulation-directory"] .. "/" .. lmrCfg["grid-directory"]
end

function lmr_config.gridMetadataFilename(id)
   -- When id is supplied give individual block
   if id then
      return lmr_config.gridDirectory() .. "/" .. lmrCfg["grid-prefix"] .. "-" .. string.format(lmrCfg["block-index-format"], id) .. lmrCfg["metadata-extension"]
   end
   -- else, return global metadataname
   return lmr_config.gridDirectory() .. "/" .. lmrCfg["grid-prefix"] .. lmrCfg["metadata-extension"]
end

function lmr_config.gridFilename(id)
   local gname = lmr_config.gridDirectory() .. "/" .. lmrCfg["grid-prefix"] .. "-" .. string.format(lmrCfg["block-index-format"], id)
   if (config.grid_format == "gziptext") then
      gname = gname .. lmrCfg["gzip-extension"]
   end
   return gname
end

function lmr_config.snapshotDirectory(snapshot)
   local dname = lmrCfg["simulation-directory"]
   dname = dname .. "/"
   dname = dname .. lmrCfg["snapshot-directory"]
   if snapshot then
      dname = dname .. "/"
      dname = dname .. string.format(lmrCfg["snapshot-index-format"], snapshot)
   end
   return dname
end

function lmr_config.fluidFilename(snapshot, blkId)
   local fname = lmr_config.snapshotDirectory(snapshot)
   fname = fname .. "/"
   fname = fname .. lmrCfg["fluid-prefix"] .. "-" .. string.format(lmrCfg["block-index-format"], blkId)
   return fname
end

function lmr_config.solidFilename(snapshot, blkId)
   local fname = lmr_config.snapshotDirectory(snapshot)
   fname = fname .. "/"
   fname = fname .. lmrCfg["solid-prefix"] .. "-" .. string.format(lmrCfg["block-index-format"], blkId)
   return fname
end

function lmr_config.gridForSimFilename(snapshot, blkId)
   local gname = lmr_config.snapshotDirectory(snapshot)
   gname = gname .. "/"
   gname = gname .. lmrCfg["grid-prefix"] .. "-" .. string.format(lmrCfg["block-index-format"], blkId)
   if (config.grid_format == "gziptext") then
      gname = gname .. lmrCfg["gzip-extension"]
   end
   return gname
end


return lmr_config

