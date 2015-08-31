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

module(..., package.seeall)

require 'lex_elems'

function transformRateConstant(t, coeffs, anonymousCollider)
   local nc = 0
   for _,c in ipairs(coeffs) do
      nc = nc + c
   end
   
   if anonymousCollider then nc = nc + 1 end

   local base = 1.0e-6
   local convFactor = base^(nc-1)

   m = {}
   m.model = t[1]

   if m.model == 'Arrhenius' then
      m.A = t.A*convFactor
      m.n = t.n
      m.C = t.C
   end

   return m
end

function rateConstantToLuaStr(rc)
   local str = ""
   if rc.model == 'Arrhenius' then
      str = string.format("{model='Arrhenius', A=%16.12e, n=%f, C=%16.12e }", rc.A, rc.n, rc.C)
   end
   return str
end

function ecModelToLuaStr(ec)
   return "{}"
end

function transformSpeciesStr(sp)
   if string.match(sp, '+') then
      return string.gsub(sp, '+', '_plus')
   end
   if sp == 'e-' then
      return 'e_minus'
   end
   -- In all other cases return string unaltered
   return sp
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
	    Mechanism = lpeg.Ct(Reaction * ( RArrow + FArrow ) * Reaction),
	    Reaction = lpeg.Ct(Participant * (Plus * Participant)^0 * (PressureDependent / pdstring)^0 ) * Space,
	    Participant = lpeg.Ct(lpeg.C(Number^0) * Space * Species * Space)
	 }

G = Space * G * -1

SpeciesG = lpeg.P{ Species }

function parseReactionString(s)
   t = lpeg.match(G, s)
   return t
end

function validateReaction(t)
   if type(t[1]) ~= 'string' then
      print("There was an error when parsing reaction number: ", #reactions+1)
      print("The first entry should be a string denoting the reaction mechanism.")
      print("Bailing out!")
      os.exit(1)
   end

   reac = parseReactionString(t[1])
   if reac == nil then
      print("There was an error parsing the reaction string for reaction number: ", #reactions+1)
      print("It seems the string is badly formed.  The given string is: ")
      print(t[1])
      print("Bailing out!")
      os.exit(1)
   end

   mass, charge = checkEquationBalances(reac, #reactions+1)
	    
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

local Species2 = lpeg.Ct(lpeg.Ct((lpeg.C(Element) * lpeg.C(Number^0) * lpeg.C(PM^0) * Solid^0))^1)
Species2G = lpeg.P{ Species2 }

function checkEquationBalances(r, rnumber)
   -- Checks both mass and charge balance
   local elem = {}
   local charge = 0
   
   for _,p in ipairs(r[1]) do
      if type(p) == 'table' then
	 local coeff = tonumber(p[1]) or 1
	 if p[2] ~= "M" then
	    -- break species into elements
	    local ets = lpeg.match(Species2, p[2])
	    for _,et in ipairs(ets) do
	       local e = et[1]
	       local ne = tonumber(et[2]) or 1
	       -- collate the element counts
	       if e ~= "e" then -- don't add electron to mass balance
                  if elem[e] then
		      elem[e] = elem[e] + coeff*ne
	          else
		      elem[e] = coeff*ne
	          end
               end
	       -- collate the charge counts
	       if et[3] == "+" then
		  charge = charge + coeff*1
	       elseif et[3] == "-" then
		  charge = charge - coeff*1
	       end
	    end
	 end
      end
   end

   -- now check on the other side of the reaction
   for _,p in ipairs(r[3]) do
      if type(p) == 'table' then
	 local coeff = tonumber(p[1]) or 1
	 if p[2] ~= "M" then
	    -- break species into elements
	    local ets = lpeg.match(Species2, p[2])
	    for _,et in ipairs(ets) do
	       local e = et[1]
	       local ne = tonumber(et[2]) or 1
	       -- collate the element counts
	       if e ~= "e" then -- don't add mass of electron to mass balance
	       	  if elem[e] then
		      elem[e] = elem[e] - coeff*ne
	          else
		      elem[e] = -coeff*ne
	          end
	       end		  
	       -- collate the charge counts
	       if et[3] == "+" then
		  charge = charge - coeff*1
	       elseif et[3] == "-" then
		  charge = charge + coeff*1
	       end
	    end
	 end
      end
   end

   -- mass check
   mass_balances = true
   for k,v in pairs(elem) do
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
--          anonymousCollider = false,
--          pressureDependent = false,
--          efficiencies = {}
-- }
--
function transformReaction(t, species, suppress_warnings)
   r = {}
   r.equation = t[1]
   r.type = "elementary"

   rs = parseReactionString(t[1])
   r.anonymousCollider = false
   -- deal with forward elements

   reactants = {}
   for _,p in ipairs(rs[1]) do
      if type(p) == 'table' then
	 coeff = tonumber(p[1]) or 1
	 p[2] = transformSpeciesStr(p[2])
	 if p[2] == "M" then
	    r.anonymousCollider = true
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
	 p[2] = transformSpeciesStr(p[2])
	 if p[2] == "M" then
	    r.anonymousCollider = true
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
      r.frc = transformRateConstant(t.fr, r.reacCoeffs, anonymousCollider)
   else
      r.frc = {model="from equilibrium constant"}
   end
   if t.br then
      r.brc = transformRateConstant(t.br, r.prodCoeffs, anonymousCollider)
   else
      r.brc = {model="from equilibrium constant"}
   end
   if t.ec then
      r.ec = t.ec
   else
      -- By default
      r.ec = {model="from thermo",iT=0}
   end

   -- Deal with efficiencies
   if r.anonymousCollider then
      r.type = "withAnonymousCollider"
      -- All efficiencies are set to 1.0
      r.efficiencies = {}
      for i=1,#species do
	 -- the { , } table needs to have 0-offset indices in it
	 r.efficiencies[i-1] = 1.0
      end
      -- Next look at the special cases
      if t.efficiencies then
	 for k,v in pairs(t.efficiencies) do
	    if not species[k] then
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
	       -- species[k] gives back a C++ index,
	       r.efficiencies[species[k]] =  v
	    end
	 end
      end
   end

   -- Look for presence of pressure dependent reaction
--[==[ NOT UPDATED : 03-Mar-2015
         TODO when pressure dependent reactions are added.
   pressure_dependent = false
   for _,p in ipairs(rs[1]) do
      if type(p) == 'string' and p == 'pressure dependent' then
	 pressure_dependent = true
      end
   end
   for _,p in ipairs(rs[3]) do
      if type(p) == 'string' and p == 'pressure dependent' then
	 pressure_dependent = true
      end
   end

   -- Deal with any efficiencies for anonymous collision partners (if needed)
   if anonymousCollider then --  or pressure_dependent then
      if anonymousCollider then
	 r.type = "withAnonymousCollider"
      end
      -- All efficiencies are set to 1.0
      r.efficiencies = {}
      for i=1,species.size do
	 -- the { , } table needs to have 0-offset indices in it
	 r.efficiencies[i] = {i-1, 1.0}
      end
      -- Next look at the special cases
      if t.efficiencies then
	 for k,v in pairs(t.efficiencies) do
	    if not species[k] then
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
	       -- species[k] gives back a C++ index,
	       -- need to insert it at +1
	       r.efficiencies[species[k]+1] = {species[k], v}
	    end
	 end
      end
      -- For efficient execution in inner loops,
      -- remove values whose efficiency is zero.
      local essentially_zero = 1.0e-9
      for i=#r.efficiencies,1,-1 do
	 if r.efficiencies[i][2] <= essentially_zero then
	    table.remove(r.efficiencies, i)
	 end
      end

      if pressure_dependent then
	 if r.frc then
	    if r.frc.model == "pressure dependent" then
	       r.frc.efficiencies = r.efficiencies
	    end
	 end
	 if r.brc then
	    if r.brc.model == "pressure dependent" then
	       r.brc.efficiencies = r.efficiencies
	    end
	 end
      end
--]==]      
   -- Look for chemistry_energy_coupling field
   --if t.chemistry_energy_coupling then
   --    r.chemistry_energy_coupling = t.chemistry_energy_coupling
   --end

   return r
end

function arrayIntsToStr(array)
   local str = "{"
   for _,val in ipairs(array) do
      str = str..string.format(" %d,", val)
   end
   str = str.."}"
   return str
end

function reacToLuaStr(r, i)
   local rstr = string.format("reaction[%d] = {\n", i)
   rstr = rstr .. string.format("  equation = \"%s\",\n", r.equation)
   rstr = rstr .. string.format("  type = \"%s\",\n", r.type)
   rstr = rstr .. string.format("  frc = %s,\n", rateConstantToLuaStr(r.frc))
   rstr = rstr .. string.format("  brc = %s,\n", rateConstantToLuaStr(r.brc))
   if r.ec then
      rstr = rstr .. string.format("  ec = %s,\n", ecModelToLuaStr(r.ec))
   end
   rstr = rstr .. string.format("  reacIdx = %s,\n", arrayIntsToStr(r.reacIdx))
   rstr = rstr .. string.format("  reacCoeffs = %s,\n", arrayIntsToStr(r.reacCoeffs))
   rstr = rstr .. string.format("  prodIdx = %s,\n", arrayIntsToStr(r.prodIdx))
   rstr = rstr .. string.format("  prodCoeffs = %s,\n", arrayIntsToStr(r.prodCoeffs))
   rstr = rstr .. string.format("  anonymousCollider = %s,\n", tostring(r.anonymousCollider))
   rstr = rstr .. "}"
   return rstr
end

function reacToStoichMatrix(spIdx, coeffs, nsp)
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

function effToStr(reac, nsp)
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
