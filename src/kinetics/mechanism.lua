-- Author: Rowan J. Gollan
-- Date: 2021-04-08
--

module(..., package.seeall)

local pow = math.pow
local sqrt = math.sqrt

require 'lex_elems'

function extractVibrationalRelaxers(mechanisms)
   local vibrational_relaxers = {}

   for p,table in pairs(mechanisms) do
      hasvib = false
      for q,mech in pairs(table) do
         if mech.type=='V-T' then hasvib = true end
      end

      if hasvib then
         vibrational_relaxers[#vibrational_relaxers+1] = p
      end
   end

   return vibrational_relaxers
end

function transformRelaxationTime(rt, p, q, db)
   local t = {}
   t.model = rt[1]
   if t.model == "Millikan-White" then
      M_p = db[p].M*1000.0 -- kg -> g
      M_q = db[q].M*1000.0 -- kg -> g
      mu = (M_p * M_q)/(M_p + M_q)
      theta_v = db[p].theta_v
      t.a = 1.16e-3*sqrt(mu)*pow(theta_v, 4/3)
      t.b = 0.015*pow(mu, 1/4)
   elseif t.model == "ParkHTC" then
      t.submodel = transformRelaxationTime(rt.submodel, p, q, db)
      if rt.sigma == nil then
         t.sigma = 1.0e-20 -- Default collision cross section in m^2 (TODO: Is this a good idea?)
      else
         t.sigma = rt.sigma
      end
   else
      print("The relaxation time model: ", t.model, " it not known.")
      print("Bailing out!")
      os.exit(1)
   end

   return t
end

function relaxationTimeToLuaStr(rt)
   local str = ""
   if rt.model == "Millikan-White" then
      str = string.format("{model='Millikan-White', a=%.3f, b=%.6f}", rt.a, rt.b)
   elseif rt.model == "ParkHTC" then
      submodelstr = relaxationTimeToLuaStr(rt.submodel)
      str = string.format("{model='ParkHTC', sigma=%.4e, submodel=%s}", rt.sigma, submodelstr)
   else
      print(string.format("ERROR: relaxation time model '%s' is not known.", rt.model))
      os.exit(1)
   end
   return str
end

function tableToString(o)
   -- Recursively evaluate a table to produce a lua-readable string
   --
   -- Notes: This function assumes that you don't have any "bad" key names
   --        such as "NO+" which would need to  be wrapped like ["NO+"]

   -- If the object o is a table, loop through it and concatenate the things inside it
   if type(o) == 'table' then
      local s = '{'

      -- We need to know how big the table is to not put a ", " at the last entry
      local tablesize = 0
      for _,__ in pairs(o) do tablesize = tablesize + 1 end

      local n = 0
      for k,v in pairs(o) do
         -- entries with a named key just get put by themselves
         if type(k) == 'number' then
            s = s .. tableToString(v)
         end
         -- Entries with a named key are written as key=value,
         if type(k) ~= 'number' then
            s = s .. k ..' = ' .. tableToString(v)
         end

         -- Possibly put a comma and a space between entries
         n = n + 1
         if n ~= tablesize then
             s = s .. ', '
         end
      end
      return s .. '}'

   -- The object o is not a table, so just return item in string form.
   else
      if type(o) == 'string' then
          return "'" .. o .. "'"
      else
          return tostring(o)
      end
   end
end

function mechanismToLuaStr(index, m)
   local typeStr
   local argStr
   if m.type == "V-T" then
      argStr = string.format("  rate = '%s',\n", m.rate)
      argStr = argStr .. string.format("  relaxation_time = %s\n", relaxationTimeToLuaStr(m.rt))
   elseif m.type == "E-T" then
      argStr = string.format("  exchange_cross_section = %s\n", tableToString(m.exchange_cross_section))
   else
      print("ERROR: type is not known: ", type)
      os.exit(1)
   end
   local mstr = string.format("mechanism[%s] = {\n", index)
   mstr = mstr .. string.format("  type = '%s',\n", m.type)
   mstr = mstr .. string.format("  p = '%s', q = '%s',\n", m.p, m.q)
   mstr = mstr .. argStr
   mstr = mstr .. "}\n"
   return mstr
end



for k,v in pairs(lex_elems) do
   _G[k] = v
end

-- Grammar

local PColliders = lpeg.V"PColliders"
local QColliders = lpeg.V"QColliders"
local Mechanism = lpeg.V"Mechanism"

G = lpeg.P{ Mechanism,
            Mechanism = lpeg.Ct( PColliders * DoubleTilde * QColliders ),
            PColliders = lpeg.Ct( Species + (Open * ( (Species * Comma^0)^1 + MolcColliders )  * Close) ),
            QColliders = lpeg.Ct( Species + (Open * ( (Species * Comma^0)^1 + AllColliders )  * Close) )
}

G = Space * G * Space * -1

function parseMechString(s)
   t = lpeg.match(G, s)
   return t
end

function validateMechanism(m, i)
   if type(m[1]) ~= 'string' then
      print("There was an error when parsing mechanism number: ", i)
      print("The first entry should be a string denoting the involved colliders.")
      print("Bailing out!")
      os.exit(1)
   end

   mech = parseMechString(m[1])
   if mech == nil then
      print("There was an error parsing the mechanism string for mechanism: ", i)
      print("It seems the string is badly formed. The given string is: ")
      print(m[1])
      print("Bailing out!")
      os.exit(1)
   end

   return true
end

function expandKeywords(keyword, species, db)
   -- Sometimes we don't want to specify every species in a list. This function defines some handy
   -- aliases for large groups of species, but will also return a table if given a valid species.
   --
   -- Examples:
   --    *molcs -> {N2, O2, NO}
   --    N2     -> {N2}
   -- Notes: If you add a new keyword to this list, also go update lex_elems.lua, and add it to the grammar bit below!
   -- @author: Nick Gibbons
   local e = {}
   if keyword=='*molcs' then
      for i,s in ipairs(species) do
         if db[s].type == 'molecule' then e[#e+1] = s end
      end

   elseif keyword=='*all' then
      for i,s in ipairs(species) do
         e[#e+1] = s
      end

   elseif keyword=='*heavy' then
      for i,s in ipairs(species) do
         if s ~= 'e-' then e[#e+1] = s end
      end

   elseif keyword=='*ions' then
      for i,s in ipairs(species) do
         if db[s].charge > 0 then e[#e+1] = s end
      end

   else
       -- If we've been given something else check that it is in the species list
      is_valid_species = false
      for i,s in ipairs(species) do
         if keyword==s then is_valid_species=true end
      end
      if is_valid_species then
           e[#e+1] = keyword
      else
         print("Error: ", keyword, "is not a valid species or expandable keyword")
         os.exit(1)
      end
   end
   return e
end

function expandColliders(t, species, db)
   -- Process the mechanism string to produce a valid table of species specified in the mechanism
   --
   -- Examples:
   --   (*molcs) -> {N2, O2, NO}
   --   (N2, O2, *molcs) -> {N2, O2, NO}
   --   (N2, O2, NO) -> {N2, O2, NO}
   --   (N2, N2, N2) -> {N2}
   -- @author: Nick Gibbons
   local ps = {}
   local unique_ps = {}

   for i,p in ipairs(t) do
      eps = expandKeywords(p, species, db)
      for j,ep in ipairs(eps) do
         if unique_ps[ep] == nil then
            unique_ps[ep] = true
            ps[#ps+1] = ep
         end
      end
   end
   return ps
end


function addUserMechToTable(index, m, mechanisms, species, db)
   t = parseMechString(m[1])
   ps = expandColliders(t[1], species, db)
   qs = expandColliders(t[2], species, db)

   for _,p in ipairs(ps) do
      for __,q in ipairs(qs) do
         mechanisms[index] = {}
         mechanisms[index].type = m.type
         mechanisms[index].p = p
         mechanisms[index].q = q
         if m.type == 'V-T' then
            mechanisms[index].rate = m.rate
            mechanisms[index].rt = transformRelaxationTime(m.relaxation_time, p, q, db)
         elseif m.type == 'E-T' then
            mechanisms[index].exchange_cross_section = m.exchange_cross_section
         end
         index = index + 1
      end
   end
   return index
end

for k,v in pairs(lex_elems) do
   _G[k] = v
end

-- Grammar

local PColliders = lpeg.V"PColliders"
local QColliders = lpeg.V"QColliders"
local Mechanism = lpeg.V"Mechanism"

G = lpeg.P{ Mechanism,
            Mechanism = lpeg.Ct( PColliders * DoubleTilde * QColliders ),
            PColliders = lpeg.Ct( Species + (Open * ( (Species * Comma^0)^1 + MolcColliders + AllColliders + HeavyColliders + IonColliders )  * Close) ),
            QColliders = lpeg.Ct( Species + (Open * ( (Species * Comma^0)^1 + MolcColliders + AllColliders + HeavyColliders + IonColliders )  * Close) )
}

G = Space * G * Space * -1









--[=[ TESTING
testStrings = {
   "N2 ~~ N2",
   "N2 ~~ (*all)",
   "O2 ~~ (N2, O2, O)",
   "(N2, O2) ~~ (*all)",
   "(*molcs) ~~ (*all)"
}

for _,tStr in ipairs(testStrings) do
   t = lpeg.match(G, tStr)
   print("p-colliders:")
   for k,v in pairs(t[1]) do print(k, v) end
   print("q-colliders:")
   for k,v in pairs(t[2]) do print(k, v) end
end
--]=]
