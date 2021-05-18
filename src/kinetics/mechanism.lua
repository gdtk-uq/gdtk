-- Author: Rowan J. Gollan
-- Date: 2021-04-08
--

module(..., package.seeall)

local pow = math.pow
local sqrt = math.sqrt

require 'lex_elems'

function transformRelaxationTime(rt, p, q, db)
   local t = {}
   t.model = rt[1]
   if t.model == "Millikan-White" then
      M_p = db[p].M*1000.0 -- kg -> g
      M_q = db[q].M*1000.0 -- kg -> g
      mu = (M_p * M_q)/(M_p + M_q)
      theta_v = db[p].theta_v
      -- If the user has not supplied values for a and b compute them from Millikan and White's correlation
      if rt.a == nil then
         t.a = 1.16e-3*sqrt(mu)*pow(theta_v, 4/3)
      else
         t.a = rt.a
      end
      if rt.b == nil then
         t.b = 0.015*pow(mu, 1/4)
      else
         t.b = rt.b
       end
   elseif t.model == "ParkHTC" then
      t.submodel = transformRelaxationTime(rt.submodel, p, q, db)
      if rt.sigma == nil then
         t.sigma = 1.0e-20 -- Default collision cross section in m^2
      else
         t.sigma = rt.sigma
      end
   elseif t.model == "ParkHTC2" then
      t.submodel = transformRelaxationTime(rt.submodel, p, q, db)
      if rt.sigma == nil then
         t.sigma = 1.0e-21 -- Default collision cross section in m^2
      else
         t.sigma = rt.sigma
      end
   elseif t.model == "KimHTC" then
      t.submodel = transformRelaxationTime(rt.submodel, p, q, db)
      t.sigma = rt.sigma
      t.exponent = rt.exponent
   else
      print("The relaxation time model: ", t.model, " it not known.")
      print("Bailing out!")
      os.exit(1)
   end

   return t
end

function calculateDissociationEnergy(dissociating_species, db)
   -- Calculate the dissociation energy by assuming the species splits
   -- completely into its atomic components. This may require some scrunity
   -- for polyatomic molecules, e.g. CO2
   -- 
   -- @author: Nick Gibbons
   local D = -db[dissociating_species].Hf
   for k,v in pairs(db[dissociating_species].atomicConstituents) do
      if db[k] == nil then
          print(string.format("ERROR: Atom %s (from %s) not present in gas model.", k, dissociating_species))
          os.exit(1)
      end
      D = D + v*db[k].Hf
   end
   return D
end

function buildChemistryCouplingModel(m, dissociating_species, db)
   -- Coupling models for chemistry vibrational coupling need their parameters
   -- assembled into a table. We do that here, computing any things that they need.
   --
   -- Notes: Consider 
   -- @author: Nick Gibbons
   local ccm = {}
   ccm.model = m.model

   if m.model == "ImpartialDissociation" then
       ccm.Thetav = db[dissociating_species].theta_v
       ccm.D = calculateDissociationEnergy(dissociating_species, db)
   else
      print(string.format("ERROR: chemistry coupling model '%s' is not known.", m.model))
      os.exit(1)
   end
   return ccm
end

function relaxationTimeToLuaStr(rt)
   local str = ""
   if rt.model == "Millikan-White" then
      str = string.format("{model='Millikan-White', a=%.3f, b=%.6f}", rt.a, rt.b)
   elseif rt.model == "ParkHTC" then
      submodelstr = relaxationTimeToLuaStr(rt.submodel)
      str = string.format("{model='ParkHTC', sigma=%.4e, submodel=%s}", rt.sigma, submodelstr)
   elseif rt.model == "ParkHTC2" then
      submodelstr = relaxationTimeToLuaStr(rt.submodel)
      str = string.format("{model='ParkHTC2', sigma=%.4e, submodel=%s}", rt.sigma, submodelstr)
   elseif rt.model == "KimHTC" then
      submodelstr = relaxationTimeToLuaStr(rt.submodel)
      str = string.format("{model='KimHTC', sigma=%.4e, exponent=%.4f, submodel=%s}", rt.sigma, rt.exponent, submodelstr)
   else
      print(string.format("ERROR: relaxation time model '%s' is not known.", rt.model))
      os.exit(1)
   end
   return str
end

function chemistryCouplingTypeToLuaStr(ct)
   local str = ""
   if ct.model == "ImpartialDissociation" then
      str = string.format("{model='ImpartialDissociation', D=%f, Thetav=%f}", ct.D, ct.Thetav)
   else
      print(string.format("ERROR: Chemistry coupling type '%s' is not known.", ct.model))
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
      argStr = string.format("  p = '%s', q = '%s',\n", m.p, m.q)
      argStr = argStr .. string.format("  rate = '%s',\n", m.rate)
      argStr = argStr .. string.format("  relaxation_time = %s\n", relaxationTimeToLuaStr(m.rt))
   elseif m.type == "E-T" then
      argStr = string.format("  p = '%s', q = '%s',\n", m.p, m.q)
      argStr = argStr .. string.format("  rate = 'ElectronExchange',\n")
      argStr = argStr .. string.format("  exchange_cross_section = %s\n", tableToString(m.exchange_cross_section))
   elseif m.type == "C-V" then
      argStr = string.format("  p = '%s',\n", m.p)
      argStr = argStr .. string.format("  rate = '%s',\n", m.rate)
      argStr = argStr .. string.format("  reaction_index = %d,\n", m.reaction_index)
      argStr = argStr .. string.format("  coupling_model = %s\n", chemistryCouplingTypeToLuaStr(m.coupling_model))
   else
      print("ERROR: mechanism type is not known: ", type)
      os.exit(1)
   end
   local mstr = string.format("mechanism[%s] = {\n", index)
   mstr = mstr .. string.format("  type = '%s',\n", m.type)
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

function seekEntriesWithLabel(table, keyname, keyvalue)
   -- Search through a table of tables, compiling any that have a key named
   -- "keyname" with a value equal to "keyvalue".
   -- 
   -- Example:
   -- t = {one = {label="reaction1", A=1, B=2}, two={label="reaction2", A=1, B=3}}
   -- out = seekEntriesWithLabel(t, "label", "reaction1"))
   -- print(tableToString(out))
   --   >> {one = {label="reaction1", A=1, B=2}}
   -- 
   -- @author: Nick Gibbons
   local entries = {}
   for k,v in pairs(table) do
      if v[keyname] == keyvalue then entries[k] = v end
   end
   return entries
end

function concatenateTables(t1, t2)
   -- Take two tables and stick them together, throwing away key information.
   -- 
   -- Example:
   -- out = concatenateTables({1,2,3}, {a=4,b=5,c=6})
   -- print(tableToString(out))
   --   >> {1,2,3,4,5,6}
   -- 
   -- @author: Nick Gibbons
   local newtable = {}
   local index = 1
   for _,v in ipairs(t1) do
      newtable[index] = v
      index = index + 1
   end
   for _,v in ipairs(t2) do
      newtable[index] = v
      index = index + 1
   end
   return newtable
end

function addUserChemMechToTable(index, m, mechanisms, species, db, reaction)
   -- Look through the surface definitions for Chemistry/Vibration coupling
   -- mechanisms. We loop through the reaction_labels table and look for
   -- reactions with matching labels, creating a deep definition of the 
   -- CV mechanism for each match. 
   --
   -- @author: Nick Gibbons
   for _,target_label in ipairs(m.reaction_labels) do

      target_reactions = seekEntriesWithLabel(reaction, "label", target_label)
      if next(target_reactions) == nil then
          print(string.format('ERROR: No reactions matching label "%s" found.', target_label))
          os.exit(1)
      end

      for reaction_index, target_reaction in pairs(target_reactions) do
         participant_species = {}
         participant_indices = concatenateTables(target_reaction.prodIdx, target_reaction.reacIdx)
         for k,v in ipairs(participant_indices) do
            spname = species[v+1] -- Lua indices start from zero, and the table is d indices
            participant_species[spname] = db[spname]
         end

         participant_molecules = seekEntriesWithLabel(participant_species, "type", "molecule")

         for speciesname,dbentry in pairs(participant_molecules) do
            mechanisms[index] = {}
            mechanisms[index].type = m.type
            mechanisms[index].p = speciesname
            mechanisms[index].rate = m.rate
            mechanisms[index].reaction_index = reaction_index - 1 -- d indexing starts at zero!
            mechanisms[index].coupling_model = buildChemistryCouplingModel(m.coupling_model, speciesname, db)
            index = index + 1
         end
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
