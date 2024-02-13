-- Module for distributing blocks to MPI tasks, required by prep.lua.
--
-- Authors: PJ and RJG
--

local mpi = {}

function mpi.mpiDistributeBlocks(args)
   -- Assign blocks to MPI tasks,
   -- keeping a record of the MPI rank (or task) for every block.
   -- This record is stored in the global variable mpiTasks.
   --
   if args and not(type(args) == "table") then
      error("mpiDistributeBlocks expects its arguments in single table with named fields", 2);
   end
   args = args or {}
   local flag = checkAllowedNames(args, {"ntasks", "dist", "preassign"})
   if not flag then
      error("Invalid name for item supplied to mpiDistributeBlocks.", 2)
   end
   --
   local nBlocks = #fluidBlocks + #solidBlocks
   -- If nTasks is not given, assume that we want one per block.
   local nTasks = (args and args.ntasks) or nBlocks
   nTasks = math.min(nTasks, nBlocks)
   --
   local option = (args and args.dist) or "round-robin"
   --
   -- The following list will eventually hold the MPI-rank for each FluidBlock.
   local mpiTaskList = {}; for i=1,nBlocks do mpiTaskList[i] = -1 end
   --
   -- Preassigned blocks.
   -- Entries in the table are of the form blockId:mpiRank
   -- Remember that blockId and mpiRank values start at zero.
   if args and (type(args.preassign) == "table") then
      for blkId,mpirank in pairs(args.preassign) do mpiTaskList[blkId+1] = mpirank end
   end
   --
   -- Distribute the rest of the blocks.
   if option == "uniform" or option == "round-robin" or option == "roundrobin" then
      -- For each block, assign to an MPI task, round-robin order.
      local mpirank = 0
      for i=1,nBlocks do
         if mpiTaskList[i] < 0 then
            -- This block not previously assigned a rank.
            mpiTaskList[i] = mpirank
            mpirank = mpirank + 1
            if mpirank == nTasks then mpirank = 0 end
         end
      end
   elseif option == "loadbalance" or option == "load-balance" then
      -- Load-balance procedure first sorts the blocks by size...
      -- first process fluidblocks
      local fluidblksNcells = {}
      local fluidtotalCells = 0
      for i, blk in ipairs(fluidBlocks) do
         fluidblksNcells[#fluidblksNcells+1] = {i, blk.ncells}
         fluidtotalCells = fluidtotalCells + blk.ncells
      end
      table.sort(fluidblksNcells, function (a,b) return a[2] > b[2] end)

      -- next process solidblocks
      local solidblksNcells = {}
      local solidtotalCells = 0
      for i, blk in ipairs(solidBlocks) do
         solidblksNcells[#solidblksNcells+1] = {i+#fluidBlocks, blk.ncells}
         solidtotalCells = solidtotalCells + blk.ncells
      end
      table.sort(solidblksNcells, function (a,b) return a[2] > b[2] end)

      -- now append both tables together
      -- we do this to ensure that every MPI task has at least one fluidblock
      -- present when performing coupled fluid-solid domain problems otherwise
      -- Eilmer balks on MPI processes that only have solidblocks present
      local blksNcells = {}
      for _, bNc in ipairs(fluidblksNcells) do blksNcells[#blksNcells+1] = bNc end
      for _, bNc in ipairs(solidblksNcells) do blksNcells[#blksNcells+1] = bNc end
      local totalCells = fluidtotalCells + solidtotalCells

      -- ...then distributes the blocks to the tasks,
      -- biggest block first into the task with the smallest load.
      -- We shall tally the loads, in number of cells, for each MPI task.
      local taskLoads = {}; for i=1,nTasks do taskLoads[i] = 0 end
      -- Account for the preassigned blocks.
      for ib=1,nBlocks do
         mpirank = mpiTaskList[ib]
         if mpirank >= 0 then
            if (ib <= #fluidBlocks) then
               taskLoads[mpirank+1] = taskLoads[mpirank+1] + fluidBlocks[ib].ncells
            else
               taskLoads[mpirank+1] = taskLoads[mpirank+1] + solidBlocks[ib - #fluidBlocks].ncells
            end
         end
      end
      -- Distribute remaining blocks.
      for _,bNc in ipairs(blksNcells) do
         local ib = bNc[1]; local ncells = bNc[2]
         if mpiTaskList[ib] < 0 then
            -- Add the so-far-unassigned block to the MPI task with smallest load.
            local indxSmallest = 1; local smallestLoad = taskLoads[1]
            for i=2,#taskLoads do
               if taskLoads[i] < smallestLoad then
                  indxSmallest = i; smallestLoad = taskLoads[i]
               end
            end
            mpiTaskList[ib] = indxSmallest-1 -- MPI task ids start from zero
            taskLoads[indxSmallest] = taskLoads[indxSmallest] + ncells
         end
      end
      --
      local maxmpiLoads = (math.max(unpack(taskLoads)))
      local minmpiLoads = (math.min(unpack(taskLoads)))
      local mpiProcessors = math.max(unpack(mpiTaskList)) + 1
      print("Load balancing - Distribute blocks to CPUs")
      print(string.format("Number of processors   \t \t = %d", mpiProcessors))
      print(string.format("Ideal cell partitioning   \t = %d cells/proc", totalCells//mpiProcessors))
      print(string.format("Smallest partition factor \t = %.3f", minmpiLoads/(totalCells/mpiProcessors)))
      print(string.format("Largest partition factor  \t = %.3f", maxmpiLoads/(totalCells/mpiProcessors)))
      print(string.format("Largest processor load    \t = %.0f cells", maxmpiLoads))

   else
      error('Did not select one of "round-robin" or "load-balance". for mpiDistributeBlocks', 2)
   end
   -- Assign the newly-constructed list to the global variable
   -- for later use in writing the job.mpimap file.
   _G.mpiTasks = mpiTaskList
   -- Finally, return the list as we have always done, however,
   -- we expect that the caller will ignore this return value.
   return mpiTaskList
end

function mpi.mpiDistributeFBArray(args)
   -- For block-marching mode, assign blocks within a FBArray object to MPI tasks,
   -- keeping a record of the MPI rank (or task) for every block.
   -- This record is stored in the global variable mpiTasks.
   --
   if args and not(type(args) == "table") then
      error("mpiDistributeFNArray expects its arguments in single table with named fields", 2);
   end
   args = args or {}
   local flag = checkAllowedNames(args, {"fba", "ntasks"})
   if not flag then
      error("Invalid name for item supplied to mpiDistributeFBArray.", 2)
   end
   if (not args.fba) or (not (args.fba.myType and (args.fba.myType == "FBArray"))) then
      error("Need to supply an FBArray object to mpiDistributeFBArray.", 2)
   end
   --
   -- If nTasks is not given, assume that we want one per block in the j,k plane.
   local nTasks = (args and args.ntasks) or (fba.njb*fba.nkb)
   -- We really want to have balanced tasks, so check that we will be allocating
   -- an equal number of blocks to each MPI task.
   if not ((args.fba.njb*args.fba.nkb)%nTasks == 0) then
      error("Number of blocks in jk plane does not divide nicely into nTasks")
   end
   -- The following list will eventually hold the MPI-rank for each FluidBlock.
   -- We may have assigned MPI tasks to other FBArray objects previously.
   local mpiTaskList = mpiTasks or {}
   -- Entries in the table are of the form blockId:mpiRank
   -- Remember that blockId and mpiRank values start at zero,
   -- but the ib, jb, kb indices start at 1 and
   -- the mpiTaskList starts at 1, to match fluidBlocks.
   taskId = 0
   if config.dimensions == 3 then
      for kb=1,args.fba.nkb do
         for jb=1,args.fba.njb do
            for ib=1,args.fba.nib do
               local blk = args.fba.blockArray[ib][jb][kb]
               mpiTaskList[blk.id+1] = taskId
            end
            taskId = taskId + 1
            if taskId >= nTasks then taskId = 0 end
         end
      end
   else
      -- Presume 2D
      for jb=1,args.fba.njb do
         for ib=1,args.fba.nib do
            local blk = args.fba.blockArray[ib][jb]
            mpiTaskList[blk.id+1] = taskId
         end
         taskId = taskId + 1
         if taskId >= nTasks then taskId = 0 end
      end
   end
   -- Assign the newly-constructed list to the global variable
   -- for later use in writing the job.mpimap file.
   _G.mpiTasks = mpiTaskList
   -- Finally, return the list as we have always done for mpiDistributeFluidBlocks,
   -- however, we expect that the caller will ignore this return value.
   return mpiTaskList
end

return mpi
