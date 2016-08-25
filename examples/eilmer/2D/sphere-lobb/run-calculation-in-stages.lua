#!/usr/bin/env dgd-lua

local floor = math.floor
coarseCellCount = 30
jobName = 'lobb'
nCells = {coarseCellCount,
          coarseCellCount*2,
	  coarseCellCount*2^2,
	  coarseCellCount*2^3 }
flowTimes = {15, 5, 3, 2}

EILMER_CMD = "e4shared"
filesToCopy = {'case.lua',
	       'air-5sp.lua',
	       'air-5sp-6r.lua',
	       'lobb.lua'}

local format = string.format

function prepare_case_file(stageIdx)
   f = assert(io.open("case.lua", "w"))
   f:write(format("ncells = %d\n", nCells[stageIdx]))
   f:write(format("no_flow_times = %f\n", flowTimes[stageIdx]))
   if stageIdx == 1 then
      f:close()
      return
   end
   -- for all other cases setup from previous case
   f:write(format("oldSoln_jobname = '%s'\n", jobName))
   f:write(format("oldSoln_dir = '../stage-%d'\n", stageIdx-1))
   -- Grab final tindx from previous simulation
   os.execute("tail -1 stage-" .. stageIdx-1 .. "/" .. jobName .. ".times > tmp")
   f1 = assert(io.open("tmp", "r"))
   final_tindx = f1:read("*number")
   f1:close()
   -- Now write to case.txt and finish
   f:write(format("oldSoln_tindx = %d\n", final_tindx))
   f:close()
   return
end

function prepare_simulation(stageIdx)
   stageDir = "stage-"..stageIdx
   -- prepare new directory
   os.execute("mkdir "..stageDir)
   -- copy over required files
   fileList = ""
   for _,f in ipairs(filesToCopy) do fileList = fileList .. f .. " " end
   os.execute("cp "..fileList.." "..stageDir.."/")
   -- run eilmer prep
   cmd = "cd "..stageDir.." && "
   cmd = cmd..EILMER_CMD.." --prep --job="..jobName
   os.execute(cmd)
   return
end

function run_simulation(stageIdx)
   stageDir = "stage-"..stageIdx
   cmd = "cd "..stageDir.." && "
   cmd = cmd..EILMER_CMD.." --run --job="..jobName
   os.execute(cmd)
   return
end

function collate_results()
   f = assert(io.open("sim-results.txt", 'w'))
   f:write("# ncells     delta     delta/D\n")
   for stageIdx, n in ipairs(nCells) do
      stageDir = "stage-"..stageIdx
      os.execute("cp shock-detachment.lua "..stageDir.."/")
      cmd = "cd "..stageDir.." && "
      cmd = cmd..EILMER_CMD.." --job="..jobName.." --custom-post --script-file=shock-detachment.lua"
      os.execute(cmd)
      f1 = io.open(stageDir.."/shock-detachment.txt"); result = f1:read("*all"); f1:close()
      f:write(string.format("%d %s", n, result))
   end
   f:close()
   return
end

function main()
   -- Run simulations in stages
   for stageIdx=1,#nCells do
      print("========================= STAGE "..stageIdx.." =============================")
      prepare_case_file(stageIdx)
      prepare_simulation(stageIdx)
      run_simulation(stageIdx)
   end
   collate_results()
   return 0
end

main()
