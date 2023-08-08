#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Remix by: Nick Gibbons
-- Date: 2021-04-08
-- History:
--    2023-02-05 Added handling for V-V rate input
--

local mechanism = require 'mechanism'

userMechs = {}
userCCMechs = {}
mechanisms = {}

local validateMechanism = mechanism.validateMechanism
local addUserMechToTable = mechanism.addUserMechToTable
local mechanismToLuaStr = mechanism.mechanismToLuaStr
local addUserChemMechToTable = mechanism.addUserChemMechToTable

function printHelp()
   print("prep-kinetics --- Prepares an energy exchange kinetics file for Eilmer.")
   print("Usage:")
   print(" > prep-kinetics gmodelFile kinInput output")
   print(" <OR> when chemistry-coupling is involved")
   print(" > prep-kinetics gmodelFile chemFile kinInput output")
   print("")
   print("   gmodelFile  : a gas model file is required as input for context")
   print("  [chemFile    : a chemistry scheme file is required as input for context when coupling models are in play]")
   print("   kinInput    : input energy exchange kinetics file in Lua format.")
   print("   output      : output file in format ready for Eilmer4.")
   print("")
   os.exit(1)
end


--%--------------------------------------------
-- Function available to user in input file.
--%--------------------------------------------

function Mechanism(t)
   -- Gather mechanisms, but don't yet validate
   userMechs[#userMechs+1] = t
end

function ChemistryCouplingMechanism(t)
   -- Gather chemistry coupling mechanisms, these don't get validated
   t.type = "C-V"
   userCCMechs[#userCCMechs+1] = t
end

-----------------------------------------------

function buildVerboseLuaFile(fName)
   f = assert(io.open(fName, 'w'))
   f:write("\n")
   f:write("mechanism = {}\n\n")
   for i,mech in ipairs(mechanisms) do
      f:write(mechanismToLuaStr(i, mech))
      f:write("\n")
   end
   f:close()
end

local function splitComponent(component)
   -- split an energy component of the form "species:energy_type", as given in the gas model
   -- into the species and energy_type
   local t = {}
   for str in string.gmatch(component, "([^:]+)") do
      table.insert(t, str)
   end
   return t[1], t[2]
end

local function buildEnergyModes(mode_names, modes)
   -- Build a table which maps the species name (and energy type) to the 
   -- energy mode index where that energy is accounted for.
   local energy_modes = {}
   for imode, mode_name in ipairs(mode_names) do
      for _, comp in ipairs(modes[mode_name]) do
         local species, energy_type = splitComponent(comp)
         if not energy_modes[species] then
            energy_modes[species] = {}
         end
         energy_modes[species][energy_type] = imode-1
      end
   end
   return energy_modes
end


function main()
   local gmodelFile, chemmodelFile, inFname, outFname
   
   if (#arg == 0 or arg[1] == "--help") then
      printHelp()
   end

   if (#arg < 3) then
      print("Not enough arguments or unknown option.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   if (#arg == 3) then
       gmodelFile = arg[1]
       inFname    = arg[2]
       outFname   = arg[3]
   end

   if (#arg == 4) then
       gmodelFile = arg[1]
       chemFile   = arg[2]
       inFname    = arg[3]
       outFname   = arg[4]
   end

   if (#arg > 4) then
      print("Too many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end


   outstring = 'Creating kinetics file "' .. outFname .. '" using input file "' .. inFname .. '"...'
   print(outstring)

   -- Execute chemical kinetics file in case we need reaction info
   -- We do this before the gas file because both have a table called
   -- "species", and we want the gas file's version to overwrite the chemFile's one.
   if chemFile ~= nil then
      dofile(chemFile)
   else
      reaction = {} 
   end

   -- Execute gas model file so we can get:
   -- 1. list of species
   -- 2. list of energy modes
   dofile(gmodelFile)

   -- The species table has indices as keys, and species names as values.
   -- Let's augment that with a reverse lookup, names as keys and indices as values.
   -- And we'll use the D-offset for species index (from 0)
   for isp,sp in ipairs(species) do
      species[sp] = isp-1
   end
   -- Do the same for energy_modes
   if energy_modes then
      energy_modes = buildEnergyModes(energy_modes, db.modes)
   elseif physical_model == "two-temperature-gas" then
      -- For 2-T, we don't require the user to set energy modes explicitly
      -- since we can make that decision. So we'll set it up.
      energy_modes = {}
      for isp,sp in ipairs(species) do
         energy_modes[sp] = {}
         energy_modes[sp].vib = 0
         energy_modes[sp].electron = 0
         energy_modes[sp].electronic = 0
      end
   elseif physical_model == "three-temperature-gas" then
      energy_modes = {}
      for isp, sp in ipairs(species) do
         energy_modes[sp] = {}
         if db[sp].type == "molecule" then
            energy_modes[sp].vib = 0
         elseif db[sp].type == "electron" then
            energy_modes[sp].electron = 1
         end
         energy_modes[sp].electronic = 1
      end
   end

   -- Load contents from user's file.
   dofile(inFname)

   index = 1
   for i,m in ipairs(userMechs) do
      if validateMechanism(m) then
         index = addUserMechToTable(index, m, mechanisms, species, db)
      else
         print("Error while trying to validate mechanism ", i)
         print(m[1])
         print("Bailing out!")
         os.exit(1)
      end
   end
   for i,m in ipairs(userCCMechs) do
      index = addUserChemMechToTable(index, m, mechanisms, species, db, reaction)
   end
   -- Now write out transformed results
   buildVerboseLuaFile(outFname)
end

main()
   


   
   
