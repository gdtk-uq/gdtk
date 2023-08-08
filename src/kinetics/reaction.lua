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
--   19-Mar-2009 :- added checking of mass balance
--                  and charge balance
--]]


local reaction = {}

local lpeg = require 'lpeg'
local lex_elems = require 'lex_elems'
local mechanism = require 'mechanism'

local calculateDissociationEnergy = mechanism.calculateDissociationEnergy

local function transformRateConstant(t, coeffs, anonymousCollider, energyModes, db)
   local nc = 0
   for _,c in ipairs(coeffs) do
      nc = nc + c
   end
   
   if anonymousCollider then nc = nc + 1 end

   local base = 1.0e-6
   local convFactor = base^(nc-1)

   local m = {}
   m.model = t[1]
   if m.model == 'Arrhenius' then
      m.A = t.A*convFactor
      m.n = t.n
      m.C = t.C
      local rctIndex = -1
      if t.rateControllingTemperature then
         if energyModes[t.rateControllingTemperature] then
            rctIndex = energyModes[t.rateControllingTemperature]
         else
            print("The supplied 'rateControllingTemperature' string is unknown: ", t.rateControllingTemperature)
            print("Bailing out!")
            os.exit(1)
         end
      end
      m.rctIndex = rctIndex
   elseif m.model == 'Arrhenius-logA' or m.model == 'Arrhenius2' then
      m.logA = t.logA
      m.B = t.B
      m.C = t.C
      local rctIndex = -1
      if t.rateControllingTemperature then
         if energyModes[t.rateControllingTemperature] then
            rctIndex = energyModes[t.rateControllingTemperature]
         else
            print("The supplied 'rateControllingTemperature' string is unknown: ", t.rateControllingTemperature)
            print("Bailing out!")
            os.exit(1)
         end
      end
      m.rctIndex = rctIndex
   elseif m.model == "fromEqConst" then
     local rctIndex = -1
     if t.rateControllingTemperature then
         if energyModes[t.rateControllingTemperature] then
            rctIndex = energyModes[t.rateControllingTemperature]
         else
            print("The supplied 'rateControllingTemperature' string is unknown: ", t.rateControllingTemperature)
            print("Bailing out!")
            os.exit(1)
         end
     end
     m.rctIndex = rctIndex
   elseif m.model == 'Park' then
      m.A = t.A*convFactor
      m.n = t.n
      m.C = t.C
      m.s = t.s
      m.mode = energyModes[t.mode] or 0
   elseif m.model == 'pressure dependent' then
      m.kInf = {}
      m.kInf.A = t.kInf.A * convFactor
      m.kInf.n = t.kInf.n
      m.kInf.C = t.kInf.C

      local convFactorLow = base^nc
      m.k0 = {}
      m.k0.A = t.k0.A * convFactorLow
      m.k0.n = t.k0.n
      m.k0.C = t.k0.C
      
      if ( t.Troe ) then
	 m.Troe = t.Troe
	 m.model = 'Troe'
      elseif ( t.YR ) then
	 m.YR = t.YR
	 m.model = 'Yungster-Rabinowitz'
      else
	 m.model = 'Lindemann-Hinshelwood'
      end
   elseif m.model == "Marrone-Treanor" then
       Runiv = 8.3145
       if t.T_D then 
          m.T_D = t.T_D
       else
          m.T_D = calculateDissociationEnergy(t.D, db) / Runiv
       end
       if t.U then
          m.U = t.U 
       elseif t.TD_on_U then
          D = calculateDissociationEnergy(t.D, db) / Runiv
          m.U = D / t.TD_on_U
       else
          print("Either U or TD_on_U needs to be provided for Marrone-Treanor reaction rate")
          print("Bailing out!")
          os.exit(1)
       end
       if t.theta_v then
          m.theta = t.theta_v
       else
          m.theta = db[t.D].theta_v or db[t.D].vib_data.theta_v
       end
       m.mode = energyModes[t.mode] or 0
       m.rate = transformRateConstant(t.rate, coeffs, anonymousCollider, energyModes, db)
   elseif m.model == "Modified-Marrone-Treanor" then
       Runiv = 8.3145
       m.T_D = t.T_D or calculateDissociationEnergy(t.D_sp, db) / Runiv
       m.theta = t.Thetav --db[t.D].theta_v
       m.aU = t.aU
       m.Ustar = t.Ustar
       m.mode = energyModes[t.mode] or 0
       m.rate = transformRateConstant(t.rate, coeffs, anonymousCollider, energyModes, db)
   else
      print("The rate constant model: ", m.model, " is not known.")
      print("Bailing out!")
      os.exit(1)
   end

   return m
end

local function rateConstantToLuaStr(rc)
   local str = ""
   if rc.model == 'Arrhenius' then
      str = string.format("{model='Arrhenius', A=%16.12e, n=%f, C=%16.12e, rctIndex=%d}", rc.A, rc.n, rc.C, rc.rctIndex)
   elseif rc.model == 'Arrhenius2' then
      str = string.format("{model='Arrhenius2', logA=%16.12e, B=%16.12e, C=%16.12e, rctIndex=%d}", rc.logA, rc.B, rc.C, rc.rctIndex)
   elseif rc.model == 'Park' then
      str = string.format("{model='Park', A=%16.12e, n=%f, C=%16.12e, s=%f, mode=%d }", rc.A, rc.n, rc.C, rc.s, rc.mode)
   elseif rc.model == 'Lindemann-Hinshelwood' then
      str = "{model='Lindemann-Hinshelwood',\n"
      str = str .. string.format("   kInf={A=%16.12e, n=%f, C=%16.12e, rctIndex=-1},\n", rc.kInf.A, rc.kInf.n, rc.kInf.C)
      str = str .. string.format("   k0={A=%16.12e, n=%f, C=%16.12e, rctIndex=-1}\n", rc.k0.A, rc.k0.n, rc.k0.C)
      str = str .. "}"
   elseif rc.model == 'Troe' then
      str = "{model='Troe',\n"
      str = str .. string.format("   kInf={A=%16.12e, n=%f, C=%16.12e, rctIndex=-1},\n", rc.kInf.A, rc.kInf.n, rc.kInf.C)
      str = str .. string.format("   k0={A=%16.12e, n=%f, C=%16.12e, rctIndex=-1},\n", rc.k0.A, rc.k0.n, rc.k0.C)
      for k,v in pairs(rc.Troe) do
	 str = str .. string.format(" %s=%16.12e, ", k, v)
      end
      str = str .. "\n}"
   elseif rc.model == 'Yungster-Rabinowitz' then
      str = "{model='Yungster-Rabinowitz',\n"
      str = str .. string.format("   kInf={A=%16.12e, n=%f, C=%16.12e, rctIndex=-1},\n", rc.kInf.A, rc.kInf.n, rc.kInf.C)
      str = str .. string.format("   k0={A=%16.12e, n=%f, C=%16.12e, rctIndex=-1},\n", rc.k0.A, rc.k0.n, rc.k0.C)
      for k,v in pairs(rc.YR) do
	 str = str .. string.format(" %s=%16.12e, ", k, v)
      end
      str = str .. "\n}"
   elseif rc.model == 'Marrone-Treanor' then
       str = "{model='Marrone-Treanor', "
       str = str .. string.format("T_D=%16.12e, U=%16.12e, theta=%16.12e, mode=%d, rate=", rc.T_D, rc.U, rc.theta, rc.mode)
       str = str .. rateConstantToLuaStr(rc.rate)
       str = str .. "}"
   elseif rc.model == 'Modified-Marrone-Treanor' then
       str = "{model='Modified-Marrone-Treanor', "
       str = str .. string.format("T_D=%16.12e, aU=%16.12e, Ustar=%16.12e, theta=%16.12e, mode=%d, rate=", rc.T_D, rc.aU, rc.Ustar, rc.theta, rc.mode)
       str = str .. rateConstantToLuaStr(rc.rate)
       str = str .. "}"
   elseif rc.model == 'fromEqConst' then
     if not rc.rctIndex then
       rc.rctIndex=-1
     end
     str = string.format("{model='fromEqConst', rctIndex=%d}", rc.rctIndex)
   else
      print(string.format("ERROR: rate constant model '%s' is not known.", rc.model))
      os.exit(1)
   end
   return str
end

local function ecModelToLuaStr(ec)
   return "{}"
end


local function checkEquationBalances(r, rnumber)
   -- Checks both mass and charge balance
   local elems = {}
   local charge = 0
   
   for _,p in ipairs(r[1]) do
      if type(p) == 'table' then
	 local coeff = tonumber(p[1]) or 1
         local sp = p[2]
	 if sp ~= "M" then
            -- look up species composition in database
            for elem,count in pairs(db[sp].atomicConstituents) do
               if elems[elem] then
                  elems[elem] = elems[elem] - coeff*count
               else
                  elems[elem] = -coeff*count
               end
	    end
            -- look up species charge in database
            charge = charge - db[sp].charge
	 end
      end
   end

   -- now check on the other side of the reaction
   for _,p in ipairs(r[3]) do
      if type(p) == 'table' then
	 local coeff = tonumber(p[1]) or 1
         local sp = p[2]
	 if p[2] ~= "M" then
            -- look up species composition in database
            for elem,count in pairs(db[sp].atomicConstituents) do
               if elems[elem] then
                  elems[elem] = elems[elem] + coeff*count
               else
                  elems[elem] = coeff*count
               end
	    end
            -- look up species charge in database
            charge = charge + db[sp].charge
	 end
      end
   end

   -- mass check
   mass_balances = true
   for k,v in pairs(elems) do
      if v ~= 0 then
	 mass_balances = false
	 print("There is a problem with the mass balance for reaction: ", rnumber)
	 print("In particular, the element: ", k)
	 print("does not balance on both sides of the reaction equation.")
      end
   end

   -- charge check
   charge_balances = true
   if charge ~= 0 then
      charge_balances = false
      print("There is a problem with the charge balance for reaction: ", rnumber)
   end

   return mass_balances, charge_balances

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
local Participant = lpeg.V"Participant"
local Reaction = lpeg.V"Reaction"
local Mechanism = lpeg.V"Mechanism"

G = lpeg.P{ Mechanism,
	    Mechanism = lpeg.Ct(Reaction * EqnSeparator * Reaction),
	    Reaction = lpeg.Ct(Participant * (Plus * Participant)^0 * (PressureDependent / pdstring)^0 ) * Space,
	    Participant = lpeg.Ct(lpeg.C(Number^0) * Space * Species * Space)
	 }

G = Space * G * -1

SpeciesG = lpeg.P{ Species }

local function parseReactionString(s)
   t = lpeg.match(G, s)
   return t
end

function reaction.validateReaction(t)
   if type(t[1]) ~= 'string' then
      print("There was an error when parsing reaction number: ", t.number)
      print("The first entry should be a string denoting the reaction mechanism.")
      print("Bailing out!")
      os.exit(1)
   end

   reac = parseReactionString(t[1])
   if reac == nil then
      print("There was an error parsing the reaction string for reaction number: ", t.number)
      print("It seems the string is badly formed.  The given string is: ")
      print(t[1])
      print("Bailing out!")
      os.exit(1)
   end

   mass, charge = checkEquationBalances(reac, t.number)
	    
   if not mass then
      print("The mass does not balance.")
      print("Bailing out!")
      os.exit(1)
   end

   if not charge then
      print("The charge does not balance.")
      print("Bailing out!")
      os.exit(1)
   end

   return true
end

-- This function transforms the user's Lua input
-- into a more verbose form for picking up in D.
--
-- An example:
-- Reaction{"H2 + I2 <=> 2 HI",
--           fr={'Arrhenius', A=1.93e14, n=0.0, C=20620.0},
--           br={'Arrhenius', A=5.4707e-7, n=0.0, C=0.0}
-- }
--
-- Gets expanded to:
-- table = {equation = "H2 + I2 <=> 2 HI",
--          type = "elementary",
--          frc = transform_rate_constant(fr),
--          brc = transform_rate_constant(fr),
--          ec = nil,
--          reacIdx = {0, 1},
--          reacCoeffs = {1, 1},
--          prodIdx = {2},
--          prodCoeffs = {2},
--          pressureDependent = false,
--          efficiencies = {}
-- }
--
function reaction.transformReaction(t, species, energyModes, db, suppressWarnings)
   r = {}
   r.equation = t[1]
   r.type = "elementary"
   r.label = t.label

   rs = parseReactionString(t[1])
   anonymousCollider = false
   -- deal with forward elements

   reactants = {}
   for _,p in ipairs(rs[1]) do
      if type(p) == 'table' then
	 coeff = tonumber(p[1]) or 1
	 if p[2] == "M" then
	    anonymousCollider = true
	 else 
	    spIdx = species[p[2]]
	    if spIdx == nil then
	       print("The following species has been declared in a reaction: ", p[2])
	       print("but is not part of the declared gas model.")
	       print("This occurred for reaction number: ", t.number)
	       if t.label then
		  print("label: ", t.label)
	       end
	       print("Bailing out!")
	       os.exit(1)
	    end
	    if reactants[spIdx] then
	       -- We've already picked this species up before
	       reactants[spIdx] = reactants[spIdx] + coeff
	    else
	       reactants[spIdx] = coeff
	    end
	 end
      end
   end
   -- We sort the table so that we get deterministic output.
   -- It helps if we need to do diffs on our D input files.
   r.reacIdx = {}
   r.reacCoeffs = {}
   local tmp = {}
   for k,v in pairs(reactants) do table.insert(tmp, k) end
   table.sort(tmp)
   for _,isp in ipairs(tmp) do
      r.reacIdx[#r.reacIdx+1] = isp
      r.reacCoeffs[#r.reacCoeffs+1] = reactants[isp]
   end

   -- do the same as above for products
   products = {}
   for _,p in ipairs(rs[3]) do
      if type(p) == 'table' then
	 coeff = tonumber(p[1]) or 1
	 if p[2] == "M" then
	    anonymousCollider = true
	 else 
	    spIdx = species[p[2]]
	    if spIdx == nil then
	       print("The following species has been declared in a reaction: ", p[2])
	       print("but is not part of the declared gas model.")
	       print("This occurred for reaction number: ", t.number)
	       if t.label then
		  print("label: ", t.label)
	       end
	       print("Bailing out!")
	       os.exit(1)
	    end
	    if products[spIdx] then
	       products[spIdx] = products[spIdx] + coeff
	    else
	       products[spIdx] = coeff
	    end
	 end
      end
   end
   r.prodIdx = {}
   r.prodCoeffs = {}
   tmp = {}
   for k,_ in pairs(products) do table.insert(tmp, k) end
   table.sort(tmp)
   for _,isp in ipairs(tmp) do
      r.prodIdx[#r.prodIdx+1] = isp
      r.prodCoeffs[#r.prodCoeffs+1] = products[isp]
   end

   if t.fr then
      r.frc = transformRateConstant(t.fr, r.reacCoeffs, anonymousCollider, energyModes, db)
   else
      r.frc = {model="fromEqConst"}
   end
   if t.br then
      r.brc = transformRateConstant(t.br, r.prodCoeffs, anonymousCollider, energyModes, db)
   else
      r.brc = {model="fromEqConst"}
   end
   if t.ec then
      r.ec = t.ec
   else
      -- By default
      r.ec = {model="from thermo",iT=0}
   end
   
   pressureDependent = false
   for _,p in ipairs(rs[1]) do
      if type(p) == 'string' and p == 'pressure dependent' then
	 pressureDependent = true
      end
   end
   for _,p in ipairs(rs[3]) do
      if type(p) == 'string' and p == 'pressure dependent' then
	 pressureDependent = true
      end
   end

   -- Deal with efficiencies
   if anonymousCollider or pressureDependent then
      if anonymousCollider then
	 r.type = "anonymous_collider"
	 r.anonymousCollider = true
      end
      -- All efficiencies are set to 1.0
      r.efficiencies = {}
      for i=0,#species-1 do
	 -- the { , } table needs to have 0-offset indices in it
	 r.efficiencies[i] = 1.0
      end
      -- Next look at the special cases
      if t.efficiencies then
	 for k,v in pairs(t.efficiencies) do
	    sp = k
	    if not species[sp] then
		  if not suppressWarnings then
		     print("WARNING: One of the species given in the efficiencies list")
		     print("is NOT one of the species in the gas model.")
		     print("The efficiency for: ", k, " will be skipped.")
		     print("This occurred for reaction number: ", t.number)
		     if t.label then
			print("label: ", t.label)
		     end
		  end
	    else
	       -- species[sp] gives back a D index,
	       r.efficiencies[species[sp]] =  v
	    end
	 end
      end
   end

   return r
end

local function arrayIntsToStr(array)
   local str = "{"
   for _,val in ipairs(array) do
      str = str..string.format(" %d,", val)
   end
   str = str.."}"
   return str
end

local function arrayNumbersToStr(array)
   local str = "{"
   for _,val in ipairs(array) do
      str = str..string.format(" %12.6e,", val)
   end
   str = str.."}"
   return str
end


function reaction.reacToLuaStr(r, i)
   local rstr = string.format("reaction[%d] = {\n", i)
   rstr = rstr .. string.format("  equation = \"%s\",\n", r.equation)
   rstr = rstr .. string.format("  type = \"%s\",\n", r.type)
   rstr = rstr .. string.format("  frc = %s,\n", rateConstantToLuaStr(r.frc))
   rstr = rstr .. string.format("  brc = %s,\n", rateConstantToLuaStr(r.brc))
   if r.ec then
      rstr = rstr .. string.format("  ec = %s,\n", ecModelToLuaStr(r.ec))
   end
   rstr = rstr .. string.format("  reacIdx = %s,\n", arrayIntsToStr(r.reacIdx))
   rstr = rstr .. string.format("  reacCoeffs = %s,\n", arrayNumbersToStr(r.reacCoeffs))
   rstr = rstr .. string.format("  prodIdx = %s,\n", arrayIntsToStr(r.prodIdx))
   rstr = rstr .. string.format("  prodCoeffs = %s,\n", arrayNumbersToStr(r.prodCoeffs))
   if r.efficiencies then
      rstr = rstr .. string.format("  efficiencies = {\n")
      -- Sort indices
      spIdx = {}
      for i,_ in pairs(r.efficiencies) do spIdx[#spIdx+1] = i end
      table.sort(spIdx)
      -- Write out efficiencies in species sorted order
      -- Skip over those efficiencies that are 0.0
      for _,i in ipairs(spIdx) do
	 if r.efficiencies[i] ~= 0.0 then
	    rstr = rstr .. string.format("    [%d]=%12.6e,\n", i, r.efficiencies[i])
	 end
      end
      rstr = rstr .. "  },\n"
   end
   if r.label then
       rstr = rstr .. string.format('  label = "%s",\n', r.label)
   end
   rstr = rstr .. "}\n"
   return rstr
end

function reaction.reacToStoichMatrix(spIdx, coeffs, nsp)
   -- Combine spIdx and coeffs into single table
   local stoichTab = {}
   for i=1,#spIdx do
      stoichTab[spIdx[i]] = coeffs[i]
   end
   local str = ""
   for isp=0,nsp-1 do
      if stoichTab[isp] then
	 str = str..string.format("%d ", stoichTab[isp])
      else
	 str = str.."0 "
      end
   end
   return str
end

function reaction.effToStr(reac, nsp)
   local str = ""
   if not reac.anonymousCollider then
      str = str.."0"
      return str
   end
   str = str.."1 "
   for isp=0,nsp-1 do
      if reac.efficiencies[isp] then
	 str = str..string.format("%12.6f ", reac.efficiencies[isp])
      else
	 str = str..string.format("%12.6f ", 0.0)
      end
   end
   return str
end

return reaction
