#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Date: 2021-04-08
--


userMechs = {}
mechanisms = {}

--%--------------------------------------------
-- Function available to user in input file.
--%--------------------------------------------

function Mechanism(t)
   -- Gather mechanisms, but don't yet validate
   userMechs[#userMechs+1]
end

-----------------------------------------------

function buildVerboseLuaFile(fName)
   f = asser(io.open(fName, 'w'))
   f:write("species = {")
   for i, sp in ipairs(species) do
      f:write(string.format("[%d]='%s', ", i-1, sp))
   end
   f:write("}\n")
   f:write("vibrational_relaxers = {")
   for i,p in mechanisms do
      f:write(string.format("%s, ", p))
   end
   f:write("}\n")
   f:write("\n")
   f:write("mechanism = {}\n")
   for _,p in ipairs(mechanisms) do
      for __,q in ipairs(mechanisms[p]) do
         f:write(mechanismToLuaStr(mechanisms[p][q], p, q))
      end
   end
   f:close()
end


function main()
   local gmodelFile, inFname, outFname
   
   if (#arg == 0 or arg[1] == "--help") then
      printHelp()
   end

   if (#arg < 3) then
      print("Not enough arguments or unknown option.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   if (#arg > 3) then
      print("Too many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   gmodelFile = arg[1]
   inFname = arg[2]
   outFname = arg[3]

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

   -- Load contents from user's file.
   dofile(inFname)

   for i,m in ipairs(userMechs) do
      if validateMechanism(m) then
         addUserMechToTable(m, mechanisms)
      else
         print("Error while trying to validate mechanism ", i)
         print(m[1])
         print("Bailing out!")
         os.exit(1)
      end
   end
   -- Now write out transformed results
   buildVerboseLuaFile(outFname)
   

end
   
   


   
   
