#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Date: 03-Mar-2009
--
-- This is mostly a port over of the chemistry
-- parser built in 2009 when I was at the NIA/NASA.
-- Here are those original notes:
--[[
-- Author: Rowan J. Gollan
-- Date: 13-Mar-2009 (Friday the 13th)
-- Place: NIA, Hampton, Virginia, USA
--
-- History:
--   20-Mar-2009 :- first put into production
--
--]]
-- In this new approach, the prep-chem program
-- is used to a take a high-level, convenient-to-use
-- input file in Lua format and convert it to a lower-level
-- format for use in the D code. Why a two-stage prep process?
-- Well, humans like to name their species with strings (eg. "N2", "O2")
-- but machines like to index into arrays with indices. Also,
-- humans like to ignore certain details of input that aren't important
-- to them, but machines don't like to guess with input.
-- Another reason is that the kinetics literature uses CGS units
-- almost exclusively, and so we allow users to input in these units.
-- The D code works in SI, so these inputs are transformed at this stage.
--

local reaction = require 'reaction'

-- We have stored some useful functions in the 'reation' module
-- Let's unpack them for use here.

local validateReaction = reaction.validateReaction
local transformReaction = reaction.transformReaction
local reacToLuaStr = reaction.reacToLuaStr
local reacToStoichMatrix = reaction.reacToStoichMatrix
local effToStr = reaction.effToStr

--%---------------------------------------------------------------------
--  Tables to store configuration 
--
--  We provide the user with some convenience function in his/her script.
--  These functions actually set entries in these tables.
--%----------------------------------------------------------------------
reactions = {}
-- Set the defaults here.
configHidden = { -- hidden from user
   tempLimits = {lower=300.0, upper=50000.0},
   odeStep = {method='rkf', errTol=1.0e-3},
   tightTempCoupling = false,
   maxSubcycles = 10000,
   maxAttempts = 4,
   __index = function (t, k) 
      return configHidden[k]
   end,
   __newindex = function (t, k, v)
      if configHidden[k] == nil then
	 print(string.format("The field '%s' cannot be set in 'config' table.", k))
      else
	 configHidden[k] = v
      end
   end,
   __call = function (_, t)
      for k, v in pairs(t) do
	 configHidden.__newindex(t, k, v)
      end
   end
}

Config = {}
setmetatable(Config, configHidden)

function odeStepToStr(o)
   if o.method == 'rkf' then
      errTol = o.errTol or 1.0e-3
      return string.format("{method='rkf', errTol=%.6e}", errTol)
   elseif o.method == 'alpha-qss' then
      eps1 = o.eps1 or 0.001
      eps2 = o.eps2 or eps1/2.0
      delta = o.delta or 1.0e-10
      maxIters = o.maxIters or 10
      str = "{method='alpha-qss', "
      str = str .. string.format("eps1= %.6e, eps2= %.6e, delta= %.6e, maxIters=%d}", eps1, eps2, delta, maxIters)
      return str
   else
      print(string.format("ERROR: The ode method name '%s' is unknown.", o.method))
      os.exit(1)
   end
end


--%----------------------------------
--  Functions available to the user
--  for configuration purposes.
--%----------------------------------

function Reaction(t)
   -- Gather reactions, but don't yet validate.
   reactions[#reactions+1] = t
end

-- This function looks for any reactions that have
-- a label which matches the supplied argument label.
-- It then removes those reactions from the list of
-- considered reactions.
function removeAllReactionsWithLabel(label)
   -- Build a list with all reactions
   -- which match label
   indices = {}
   for i,r in ipairs(reactions) do
      if label == 'string' then
	 if r.label == label then
	    indices[#indices+1] = i
	 end
      else
	 for _,l in ipairs(label) do
	    if r.label == l then
	       indices[#indices+1] = i
	    end
	 end
      end
   end

   table.sort(indices, function(a,b) return a > b end)
   for _,i in ipairs(indices) do
      table.remove(reactions, i)
   end
end

-- This function removes a reactions from the list based on
-- a number or list of numbers.
-- The global table of reactions is altered.
-- Be careful with repeated invocations of this function
-- as the positions of reactions in the table will change.
-- You may not be removing what you think you are removing
-- if you don't consider the change in table position after
-- repeated calls to this function.
function removeReaction(n)
   if type(n) == 'number' then
      table.remove(reactions, n)
   else
      table.sort(n, function(a,b) return a > b end)
      for _,i in ipairs(n) do
	 table.remove(reactions, i)
      end
   end
end

-- This is an exclusive selection. It will select only
-- the reactions with a label that matches either the
-- supplied string or strings in the list.
function selectOnlyReactionsWithLabel(t)
   -- 1. Make a list of ALL available labels
   excluded_labels = {}
   for _,r in ipairs(reactions) do
      excluded_labels[#excluded_labels+1] = r.label or 'no label'
   end
   -- 2. Build a list of those labels we are keeping
   --    that is, selecting.
   --    We will remove the 'selected' ones from the
   --    the list of ALL labels
   indices_to_remove = {}
   if type(t) == 'string' then
      for i,e in ipairs(excluded_labels) do
	 if e == t then
	    indices_to_remove[#indices_to_remove+1] = i
	 end
      end
   else
      for _,l in ipairs(t) do
	 for i,e in ipairs(excluded_labels) do
	    if e == l then
	       indices_to_remove[#indices_to_remove+1] = i
	    end
	 end
      end
   end
   -- 3. Now remove the selected ones from the list
   --    What remains is what we are removing
   table.sort(indices_to_remove, function (a,b) return a > b end)

   for _,i in ipairs(indices_to_remove) do
      table.remove(excluded_labels, i)
   end
   -- 4. Now pass off the non-selected ones for removal
   removeAllReactionsWithLabel(excluded_labels)
end

-- This exclusive selection works with a number
-- or list of numbers.
function selectOnlyReactions(n)
   -- If n is a single integer, make into a table
   -- with one entry
   if type(n) == 'number' then
      n = {n}
   end

   excluded_reactions = {}
   for i=1,#reactions do
      excluded_reactions[#excluded_reactions+1] = i
   end

   table.sort(n, function (a,b) return a > b end)

   for _,i in ipairs(n) do
      table.remove(excluded_reactions, i)
   end
   -- Now pass of all those not selected to the
   -- removal function
   removeReactions(excluded_reactions)
end

function printHelp()
   print("prep-chem --- Prepares a chemistry input file for Eilmer4.")
   print("Usage:")
   print(" > prep-chem [--compact] gmodelfile cheminput output")
   print("")
   print("   gmodelfile  : a gas model file is required as input for context")
   print("   cheminput   : input chemistry file in Lua format.")
   print("   output      : output file in format ready for Eilmer4.")
   print("")
   print("Options:")
   print("   --compact   : produce a text file called 'chem-compact-notation.inp'")
   print("                 which is used to configure a GPU chemistry kernel.")
   os.exit(1)
end

function buildVerboseLuaFile(fName)
   f = assert(io.open(fName, 'w'))
   -- Write out species in order so that we may later check for consistency with gas model
   f:write("species = {")
   for i, sp in ipairs(species) do
      f:write(string.format("[%d]='%s', ", i-1, sp))
   end
   f:write("}\n")
   -- Write out configuration settings
   f:write(string.format("config = {\n"))
   f:write(string.format("  tempLimits = {lower=%f, upper=%f},\n",
			 configHidden.tempLimits.lower,
			 configHidden.tempLimits.upper))
   f:write(string.format("  odeStep = %s,\n", odeStepToStr(configHidden.odeStep)))
   f:write(string.format("  tightTempCoupling = %s,\n", tostring(configHidden.tightTempCoupling)))
   f:write(string.format("  maxSubcycles = %d,\n", configHidden.maxSubcycles))
   f:write(string.format("  maxAttempts = %d\n", configHidden.maxAttempts))
   f:write("}\n")
   f:write("\n")
   f:write("reaction = {}\n")
   for i,r in ipairs(reactions) do
      f:write(reacToLuaStr(r, i))
      f:write("\n")
   end
   f:close()
end

function buildCompactNotationFile(fName)
   f = assert(io.open(fName, 'w'))
   f:write(string.format("%d %d\n", #species, #reactions))
   f:write("#Reactant Stoichiometric Matrix\n")
   for _,r in ipairs(reactions) do
      f:write(reacToStoichMatrix(r.reacIdx, r.reacCoeffs, #species))
      f:write("\n")
   end
   f:write("#Product Stoichiometric Matrix\n")
   for _,r in ipairs(reactions) do
      f:write(reacToStoichMatrix(r.prodIdx, r.prodCoeffs, #species))
      f:write("\n")
   end
   f:write("#Anonymous collider? No : 0, Yes : 1\n")
   for _,r in ipairs(reactions) do
      f:write(effToStr(r, #species))
      f:write("\n")
   end
   f:close()
end

--%----------------------------------
--  Main coordinating routine.
--%----------------------------------
function main()
   local gmodelName, inFname, outFname
   local doCompact = false
   if ( #arg == 0 or arg[1] == "--help" ) then
      printHelp()
   end
   
   if ( #arg < 3 ) then
      print("Not enough arguments or unknown option.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   if ( #arg > 4 ) then
      print("Too many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   if ( #arg == 4 ) then
      -- Check we did ask for compact notation.
      if ( arg[1] ~= "--compact" ) then
	 print(format.string("The option '%s' is not recognised."))
	 print("Exiting program without doing anything.")
	 printHelp()
      end
      gmodelName = arg[2]
      inFname = arg[3]
      outFname = arg[4]
      doCompact = true
   else
      gmodelName = arg[1]
      inFname = arg[2]
      outFname = arg[3]
   end
   -- Execute gas model file, just so we can get the list of species and energy modes
   dofile(gmodelName)
   -- Now we'll make the species table give us reverse look up.
   -- For a given species name, we want its D (0-offset) index
   for i, sp in ipairs(species) do
      species[sp] = i-1
   end

   if not energyModes then
      -- I didn't realise some gas models already had energyModes
      -- when I added energy_modes to the multi-T gas model...
      -- It would be nice to consolidate these in the future - RW 22/3/23
      if energy_modes then
	 energyModes = energy_modes
      else
	 energyModes = {}
      end
   end

   for i, mode in ipairs(energyModes) do
      energyModes[mode] = i-1
   end
   energyModes["translational"] = -1

   -- For two-temperature model, "db" might be empty.
   -- In which case, we'll load the entire database so that we have access
   -- to the atomicConsituents field for each species.
   if not db then
      DGD = os.getenv("DGD")
      dir = DGD.."/data/"
      dbName = dir.."species-database.lua"
      dofile(dbName)
   end

   -- Transform reactions internally
   dofile(inFname)
   for i,r in ipairs(reactions) do
      r.number = i
      if validateReaction(reactions[i]) then 
         reactions[i] = transformReaction(r, species, energyModes, db, SUPPRESS_WARNINGS)
      else
         print("Error while trying to validate reaction ", i)
         print(reactions[i][1])
         print("Bailing out!")
         os.exit(1)
      end
   end
   -- Now write out transformed results
   buildVerboseLuaFile(outFname)
   if doCompact then
      buildCompactNotationFile('chem-compact-notation.inp')
   end
end

main()
