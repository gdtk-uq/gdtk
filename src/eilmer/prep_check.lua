-- Module for checking and setting some configuration, used by prep.lua.
--
-- Authors: PJ, RJG, Kyle, Nick, Lachlan and Daryl
--
-- 2021-04-17: extracted from prep.lua

local prep_check = {}


function prep_check.initTurbulence(fs, turbulence_model_name)
    -- Setup dynamic turbulent primitive array and check user inputs
    if turbulence_model_name == "none" then
        -- Check to ensure the user hasn't tried defining any turbulence stuff
        if fs.tke~=0.0 then error(string.format("Turbulence model is none but tke set: %.18e", fs.tke)) end
        if fs.omega~=1.0 then error(string.format("Turbulence model is none but omega set: %.18e", fs.omega)) end
        if fs.nuhat~=0.0 then error(string.format("Turbulence model is none but nuhat set: %.18e", fs.nuhat)) end
        turb = {}
    elseif string.find(turbulence_model_name , "k_omega") then
        if fs.nuhat~=0.0 then error(string.format("Turbulence model is k_omega but nuhat set: %.18e", fs.nuhat)) end
        turb = {fs.tke, fs.omega}
    elseif string.find(turbulence_model_name , "spalart_allmaras") then
        if fs.tke~=0.0 then error(string.format("Turbulence model is spalart_allmaras but tke set: %.18e", fs.tke)) end
        if fs.omega~=1.0 then error(string.format("Turbulence model is spalart_allmaras but omega set: %.18e", fs.omega)) end
        turb = {fs.nuhat}
    else
        error(string.format("Unsupported turbulence model: %s", turbulence_model_name))
    end
    return turb
end

function prep_check.checkCellVolumes(t)
   if not t then
      t = {}
   end
   -- Stop reporting cells after this limit
   local badCellLimit = t.badCellLimit or 20
   local badCells = {}
   local badCellCount = 0
   for ib,blk in ipairs(fluidBlocks) do
      local grid = blk.grid
      for idx=0,grid:get_ncells()-1 do
	 local vol = grid:cellVolume(idx)
	 if vol <= 0 then
	    badCellCount = badCellCount + 1
	    badCells[#badCells+1] = {ib, idx}
	    if badCellCount >= badCellLimit then
	       return false, badCells
	    end
	 end
      end
   end
   if #badCells > 0 then
      return false, badCells
   else
      return true, badCells
   end
end

function prep_check.perform_spatial_gradient_consistency_check()
   -- Not all spatial gradient options are available, depending on the type of grid.
   -- First, search for any unstructured grids, since these are the most restricted.
   unstructuredGridsPresent = false
   for _,blk in ipairs(fluidBlocks) do
      -- in prep-flow.lua, we may not have a real grid with in the FluidBlock
      if blk.grid and (blk.grid:get_type() == "unstructured_grid") then
	 unstructuredGridsPresent = true
	 break
      end
   end
   -- In prep.lua, we don't actually make proper use of the Grid objects
   -- so guard the following tests.
   if gridsList and #gridsList > 1 then
      for _,gridMetadata in ipairs(gridsList) do
         if gridMetadata.dimensions and config.dimensions ~= gridMetadata.dimensions then
            local msg = string.format("Mismatch in dimensions, config %d grid %d.",
                                      config.dimensions, gridMetadata.dimensions)
            error(msg)
         end
      end
   end
   if unstructuredGridsPresent then
      if config.spatial_deriv_calc == "divergence" then
	 print("NOTE: config.spatial_deriv_calc is being set to 'least_squares' because unstructured grids detected.")
	 config.spatial_deriv_calc = "least_squares"
      end
      if config.spatial_deriv_locn == "vertices" then
	 print("NOTE: config.spatial_deriv_locn is being set to 'cells' when using least squares.")
	 config.spatial_deriv_locn = "cells"
      end
   else
      -- Only structured grids are present.
      -- 2D structured grids have all options available.
      if config.dimensions == 3 then
         if config.spatial_deriv_calc == "divergence" then
            print("NOTE: config.spatial_deriv_calc is being set to 'least_squares' for 3D simulations.")
            config.spatial_deriv_calc = "least_squares"
         end
      end
   end
end

function prep_check.warn_if_blocks_not_connected()
   -- It would be very unusual to have defined multiple FluidBlocks
   -- and not have any connections between them.
   -- Such an arrangement would be unintended, almost certainly.
   if #fluidBlocks > 1 then
      local n = 0
      for _,blk in ipairs(fluidBlocks) do
         for bndry,bc in pairs(blk.bcList) do
            if string.find(bc.type, "exchange_") then
               n = n + 1
            end
         end
      end
      if n == 0 then
         print("WARNING: Did not find any inter-block connections.")
         print("         Did you forget to call identifyBlockConnections()?")
      end
   end
end

function prep_check.check_DFT_settings()
    -- Check to see that the DFT has been correctly configured.

    if not ((config.DFT_n_modes * config.DFT_step_interval) == config.max_step) then
        print("WARNING: config.DFT_n_modes * config.DFT_step_interval should equal config.max_step")
    end
    if not config.fixed_time_step then
        print("WARNING: should turn config.fixed_time_step on when computing the DFT")
    end
end

return prep_check
