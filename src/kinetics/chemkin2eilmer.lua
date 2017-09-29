#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Date: 2017-09-23
--
-- This program is used to convert a chemistry scheme 
-- in the CHEMKIN input format into the chemistry input
-- format used by Eilmer.
--
-- Usage:
-- > chemkin2eilmer file-in-chemkin-fmt.inp chem-scheme.lua
-- where:
--    file-in-chemking-fmt.inp : a text file, in the chemkin format
--    chem-scheme.lua : the output Lua file for use in eilmer
--

require 'lex_elems'
require 'reaction'

local R_U_Cal = 1.9872

function split (s, sep)
  sep = lpeg.P(sep)
  local elem = lpeg.C((1 - sep)^0)
  local p = lpeg.Ct(elem * (sep * elem)^0)   -- make a table capture
  return lpeg.match(p, s)
end

-- lexical elements for parsing the whole reaction string
-- get common elements from lex_elems.lua
for k,v in pairs(lex_elems) do
   _G[k] = v
end
-- module-specific elements
local PressureDependent = Open * Plus * "M" * Space * Close
local function pdstring() return "pressure dependent" end

-- Grammar
local Comment = lpeg.P("!")
local Participant = lpeg.V"Participant"
local Reaction = lpeg.V"Reaction"
local Mechanism = lpeg.V"Mechanism"
local ForwardRate = lpeg.V"ForwardRate"
local MechanismWithForwardRate = lpeg.V"MechanismWithForwardRate"


R1 = lpeg.P{MechanismWithForwardRate,
	    MechanismWithForwardRate = lpeg.Ct(lpeg.C(Mechanism) * Space * ForwardRate * Space);
	    Mechanism = Reaction * Space * ( FRArrow + FArrow ) * Space * Reaction;
	    Reaction = Participant * (Plus * Participant)^0 * (PressureDependent)^0;
	    Participant = Number^0 * Space * Species * Space;
	    ForwardRate = lpeg.Ct(lpeg.C(Number) * Space * lpeg.C(Number) * Space * lpeg.C(Number));
	   }

R1 = Space * R1 * -1

function parseMechanismWithRate(str)
   t = lpeg.match(R1, str)
   return t
end

function printHelp()
   print("chemkin2eilmer -- Converts a Chemkin format input file into gas and chemistry files ready for eilmer.")
   print("Usage:")
   print("> chemkin2eilmer chemkin.inp gas-model.lua chem.lua")
   print("")
   print("   chemkin.inp   : a Chemkin format input file.")
   print("   gas-model.lua : (output) name for Eilmer gas model file.")
   print("   chem.lua      : (output) name for Eilmer chemistry files.")

end

local function split_string(str)
   tokens = {}
   for tk in string.gmatch(str, "%S+") do
      tokens[#tokens+1] = tk
   end
   return tokens
end

function parseChemkinFileForElements(f)
   inElementsSection = false
   while inElementsSection == false do
      line = f:read("*line")
      if not line then
	 print("End of file encountered before finding 'ELEMENTS' section.")
	 print("Exiting.")
	 os.exit(1)
      end
      tks = split_string(line)
      if tks[1] == "ELEMENTS" then
	 inElementsSection = true
      end
   end

   while inElementsSection do
      line = f:read("*line")
      if not line then
	 print("End of file encountered before reaching end of 'ELEMENTS' section.")
	 print("Exiting.")
	 os.exit(1)
      end
      -- Not really interested in picking up elements at this
      -- stage. Eilmer just needs the species.
      tks = split_string(line)
      if tks[1] == "END" then
	 inElementsSection = false
	 break
      end
   end
   return
end

function parseChemkinFileForSpecies(f)
   species = {}
   inSpeciesSection = false

   while inSpeciesSection == false do
      line = f:read("*line")
      if not line then
	 print("End of file encountered before finding 'SPECIES' section.")
	 print("Exiting.")
	 os.exit(1)
      end
      tks = split_string(line)
      if #tks > 0 then
	 if tks[1]:sub(1,4) == 'SPEC' then
	    inSpeciesSection = true
	 end
	 if #tks > 1 then
	    table.remove(tks, 1)
	    -- For all other cases, we should have legitimate species as tokens
	    for _,tk in ipairs(tks) do
	       species[#species+1] = lpeg.match(lpeg.C(Species), tk)
	    end
	 end
      end
   end

   while inSpeciesSection do
      line = f:read("*line")
      if not line then
	 print("End of file encountered for reaching end of 'SPECIES' section.")
	 print("Exiting.")
	 os.exit(1)
      end
      tks = split_string(line)
      if #tks > 0 then
	 if tks[1] == 'END' then
	    inSpeciesSection = false
	    break
	 end
	 for _,tk in ipairs(tks) do
	    if tk == 'END' then
	       inSpeciesSection = false
	       break
	    end
	    -- For all other cases, we should have legitimate species as tokens
	    sp = lpeg.match(lpeg.C(Species), tk)
	    if not sp then
	       print("Error trying to match species in 'SPECIES' section.")
	       print("Bad string: ", tk)
	       error("Syntax error", 2)
	    end
	    species[#species+1] = lpeg.match(lpeg.C(Species), tk)
	 end
      end
   end
   return species
end

function parseChemkinFileForReactions(f)
   reactions = {}
   inReactionsSection = false
   while inReactionsSection == false do
      line = f:read("*line")
      if not line then
	 print("End of file encountered before finding 'REACTIONS' section.")
	 print("Exiting.")
      end
      tks = split_string(line)
      if tks[1] == "REACTIONS" then
	 inReactionsSection = true
      end
   end

   while inReactionsSection do
      line = f:read("*line")
      if not line then
	 print("End of file encountered before reaching end of 'REACTIONCS' section.")
	 print("Exiting.")
	 os.exit(1)
      end
      print("line= ", line)
      tks = split_string(line)
      if tks[1] == 'END' then
	 inReactionsSection = false
	 break
      end

      -- Try to determine if we have a:
      -- 1. a reaction line
      -- 2. a continuation line associated with the reaction before
      -- 3. a comment line

      -- Break apart line and throw away anything with a comment, then reassemble.
      iC = -1
      for i,tk in ipairs(tks) do
	 if tk:sub(1,1) == "!" then
	    iC = i -- index of comment start
	    break
	 end
      end
      if iC > 0 then
	 line = ""
	 for i=1,iC-1 do
	    line = line .. tks[i] .. "  "
	 end
      end
      
      -- Attempt to parse a reaction string
      reac = parseMechanismWithRate(line)
      if reac then
	 reactions[#reactions+1] = {}
	 reactions[#reactions].mechanism = reac[1]
	 reactions[#reactions].forwardRate = reac[#reac]
      else
	 -- We really should find a line with extra info
	 -- about the earlier reaction line
	 reacInfo = split(line, "/")
	 if reacInfo then
	    reactions[#reactions].extraInfo = reacInfo
	 end
      end
   end
   return reactions
end

function transformReactions(reactions)
   for i,reac in ipairs(reactions) do
      reac[1] = reac.mechanism
      FR = reac.forwardRate
      print("FR[1]= ", FR[1])
      print("tonumber= ", tonumber(FR[1]))
      reac.fr = {'Arrhenius', A=tonumber(FR[1]), n=tonumber(FR[2]), C=tonumber(FR[3])/R_U_Cal}
   end
end

function writeReactionsToFile(fname, reactions)
   of = assert(io.open(fname, 'w'))
   of:write(string.format("-- Auto-generated by chemkin2eilmer on: %s\n\n",
			 os.date("%d-%b-%Y %X")))
   of:write("\n")
   for i,reac in ipairs(reactions) do
      of:write("Reaction{\n")
      of:write(string.format("        '%s',\n", reac.mechanism))
      of:write(string.format("        fr={'%s', A=%12.6e, n=%12.6e, C=%12.6e},\n",
			     reac.fr[1], reac.fr.A, reac.fr.n, reac.fr.C))
      of:write("}\n")
      of:write("\n")
   end
   of:close()
end

function main()
   if (#arg == 0 or arg[1] == "--help") then
      printHelp()
   end

   if (#arg ~= 3) then
      print("Three arguments expected.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   inFname = arg[1]
   outGasFile = arg[2]
   outChemFile = arg[3]

   -- -----------------------------------------
   -- 1. Create gas model file
   -- -----------------------------------------
   f = assert(io.open(inFname, 'r'))
   parseChemkinFileForElements(f)
   species = parseChemkinFileForSpecies(f)
   -- Try to create a gas model before moving on
   tmpFile = assert(io.open('gas-tmp.inp', 'w'))
   tmpFile:write("model = 'thermally perfect gas'\n")
   speciesStr = "species = { "
   for _,sp in ipairs(species) do
      speciesStr = speciesStr .. string.format("'%s', ", sp)
   end
   speciesStr = speciesStr .. "}\n"
   tmpFile:write(speciesStr)
   tmpFile:close()
   cmd = string.format("prep-gas gas-tmp.inp %s", outGasFile)
   os.execute(cmd)
   
   -- -------------------------------------------
   -- 2. Create reactions file
   -- -------------------------------------------
   reactions = parseChemkinFileForReactions(f)
   transformReactions(reactions)
   writeReactionsToFile(outChemFile, reactions)
   
   

   
   



end

main()
