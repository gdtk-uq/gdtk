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

local lex_elems = require 'lex_elems'
local reaction = require 'reaction'
local lpeg = require 'lpeg'

local R_U_Cal = 1.9872

function split (s, sep)
  sep = lpeg.P(sep)
  local elem = lpeg.C((1 - sep)^0)
  local p = lpeg.Ct(elem * (sep * elem)^0)   -- make a table capture
  return lpeg.match(p, s)
end
function trim(s)
  -- from PiL2 20.4
  return (s:gsub("^%s*(.-)%s*$", "%1"))
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
	    Mechanism = Reaction * Space * ( FRArrow + FArrow + Equals ) * Space * Reaction;
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
      local line = f:read("*line")
      if not line then
	 print("End of file encountered before reaching end of 'REACTIONS' section.")
	 print("Exiting.")
	 os.exit(1)
      end
      tks = split_string(line)
      if tks[1] == 'END' then
	 inReactionsSection = false
	 break
      end

      -- Try to determine if we have a:
      -- 1. a reaction line
      -- 2. a continuation line associated with the reaction before
      -- 3. a line marking a DUPLICATE reaction
      -- 4. a comment line

      -- Look for complete comment line first and ignore, so only proceed if not comment line.
      if line:sub(1,1) ~= "!" then
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
            reactions[#reactions].mechanism = trim(reac[1])
            reactions[#reactions].forwardRate = reac[#reac]
            reactions[#reactions].extraInfo = {}
         else
            -- We really should find a line with extra info
            -- about the earlier reaction line
            reacInfo = split(line, "/")

            if type(reacInfo) == 'table' then
               reacInfo[#reacInfo] = nil
               -- Remove any tokens after a comment is found.
               iC = -1
               for i,tk in ipairs(reacInfo) do
                  if tk:sub(1,1) == "!" then
                     iC = i
                     break
                  end
               end
               if iC > 0 then
                  for i=iC,#reacInfo do
                     table.remove(reacInfo)
                  end
               end

               for _,entry in ipairs(reacInfo) do
                  reactions[#reactions].extraInfo[#(reactions[#reactions].extraInfo)+1] = entry
               end

            end
         end
      end
   end
   return reactions
end

function transformReactions(reactions)
   for _,reac in ipairs(reactions) do
      reac.efficiencies = {}
      reac[1] = reac.mechanism
      FR = reac.forwardRate
      reac.fr = {'Arrhenius', A=tonumber(FR[1]), n=tonumber(FR[2]), C=tonumber(FR[3])/R_U_Cal}
      -- Transform the extra info
      -- The extra info comes in pairs of entries (1,2), (3,4), etc.
      for i=1,#reac.extraInfo,2 do
	 key = trim(reac.extraInfo[i])
	 val = trim(reac.extraInfo[i+1])
	 species = lpeg.match(Species, key)
	 -- Note: the pattern matching will actually match against 'REV' as
	 -- a legitimate species name, it thinks it's a species composed of
	 -- elements 'R', 'E', and 'V'. That's why we test and action on 
	 -- 'REV' first before trying to use the token as a species.
	 if key == 'REV' then
	    tks = split_string(val)
	    reac.br = {'Arrhenius', A=tonumber(tks[1]), n=tonumber(tks[2]), C=tonumber(tks[3])/R_U_Cal}
	 elseif key == 'LOW' then
	    tks = split_string(val)
	    reac.fr = {'pressure dependent',
		       kInf={A=reac.fr.A, n=reac.fr.n, C=reac.fr.C},
		       k0={A=tonumber(tks[1]), n=tonumber(tks[2]), C=tonumber(tks[3])/R_U_Cal}
	    }
	 elseif key == 'TROE' then
	    tks = split_string(val)
	    reac.fr.Troe = {a=tonumber(tks[1]), T3=tonumber(tks[2]), T1=tonumber(tks[3])}
	    if tks[4] then
	       reac.fr.Troe.T2 = tonumber(tks[4])
	    end
	 elseif species then
	    reac.efficiencies[#(reac.efficiencies)+1] = {species, tonumber(val)}
	 end
      end
   end
end

function writeReactionsToFile(fname, reactions)
   of = assert(io.open(fname, 'w'))
   of:write(string.format("-- Auto-generated by chemkin2eilmer on: %s\n\n",
			 os.date("%d-%b-%Y %X")))
   of:write("Config{\n")
   of:write("   odeStep = {method='alpha-qss'}\n")
   of:write("}\n")

   for i,reac in ipairs(reactions) do
      of:write("Reaction{\n")
      of:write(string.format("        '%s',\n", reac.mechanism))
      if reac.fr[1] == 'Arrhenius' then
	 of:write(string.format("        fr={'%s', A=%12.6e, n=%12.6e, C=%12.6e},\n",
				reac.fr[1], reac.fr.A, reac.fr.n, reac.fr.C))
      else
	 -- Assume pressure dependent
	 of:write(string.format("        fr={'%s',\n", reac.fr[1]))
	 of:write(string.format("            kInf={A=%12.6e, n=%12.6e, C=%12.6e},\n",
				reac.fr.kInf.A, reac.fr.kInf.n, reac.fr.kInf.C))
	 of:write(string.format("            k0={A=%12.6e, n=%12.6e, C=%12.6e},\n",
				reac.fr.k0.A, reac.fr.k0.n, reac.fr.k0.C))
	 if reac.fr.Troe then
	    of:write(string.format("            Troe={a=%12.6e, T3=%12.6e, T1=%12.6e, ",
				   reac.fr.Troe.a, reac.fr.Troe.T3, reac.fr.Troe.T1))
	    if reac.fr.Troe.T2 then
	       of:write(string.format("T2=%12.6e ", reac.fr.Troe.T2))
	    end
	    of:write("},\n")
	 end
	 of:write("        },\n")
      end
      if reac.br then
	 of:write(string.format("        br={'%s', A=%12.6e, n=%12.6e, C=%12.6e},\n",
				reac.br[1], reac.br.A, reac.br.n, reac.br.C))
      end
      if #(reac.efficiencies) > 0 then
	 of:write("        efficiencies={")
	 for _,eff in ipairs(reac.efficiencies) do
	    of:write(string.format("%s=%12.6e, ", eff[1], eff[2]))
	 end
	 of:write("},\n")
      end
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
   -- 0. Scrub input chemkin file
   -- -----------------------------------------
   -- There are some idiosyncracies allowed in Chemkin input
   -- that are troublesome to parse. We use sed to clean up
   -- the input before starting to process the file.
   --
   -- Specifically, we will:
   -- 1. Convert "AR" as Argon symbol to the proper "Ar".
   -- 2. We will NOT allow participants in a reaction to 
   --    crowd around the plus sign. For example:
   --    N2+Ar  should transform to N2 + Ar
   --    The former causes problem as the lpeg parser would
   --    try to interpret that as N2+ species but then find
   --    that it's a malformed species symbol.
   -- 3. Similarly, we'll add space around third bodies
   --    that designate a pressure-dependent reaction.
   --    OH(+M) should transform to OH (+M)
   -- 4. Convert species with (S) to _S
   --
   -- We handle these with an invocation of 'sed' and
   -- put the result in a temporary file, tmp.
   sedCmd = string.format('sed -e "s/AR/Ar/g" -e "s/\\([A-Z0-9]\\)+\\([A-Z]\\)/\\1 + \\2/g" -e "s/(+M)/ (+M)/g" -e "s/(S)/_S/g" < %s > tmp', inFname)
   os.execute(sedCmd)

   -- -----------------------------------------
   -- 1. Create gas model file
   -- -----------------------------------------
   print("Creating gas model file: ", outGasFile)
   f = assert(io.open('tmp', 'r'))
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
   tmpFile:write("options = {database='prefer-grimech'}\n")
   tmpFile:close()
   cmd = string.format("prep-gas gas-tmp.inp %s", outGasFile)
   os.execute(cmd)
   os.execute("rm gas-tmp.inp")
   print("Done: gas model file created.")
   
   -- -------------------------------------------
   -- 2. Create reactions file
   -- -------------------------------------------
   print("")
   print("Creating chemistry file for input to prep-chem: ", outChemFile)
   reactions = parseChemkinFileForReactions(f)
   transformReactions(reactions)
   writeReactionsToFile(outChemFile, reactions)
   print("Done: chemistry file created.")
   
   -- Clean up tmp file
   os.execute("rm tmp")

end

main()
