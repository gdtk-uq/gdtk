-- Author: Rowan J. Gollan
-- Date: 2021-04-08
-- 

local mechanism = {}

local pow = math.pow
local sqrt = math.sqrt

local lpeg = require 'lpeg'
local lex_elems = require 'lex_elems'

local function SSHSigma(p, db)
   -- use 4.2 as the sigma value for N2 to better match experiment
   -- as per Thivet et al.
   if p == "N2" then
      return 4.2
   else
      return db[p].sigma
   end
end

local function computeMassFactor(p, db)
   local m = {}
   for k, v in pairs(db[p].atomicConstituents) do
      for _=1,v do
         m[#m+1] = db[k].M
      end
   end
   Cp2 = (m[1]*m[1] + m[2]*m[2]) / (2*m[1]*m[2]*(m[1] + m[2]))
   return 1. / (db[p].M * Cp2)
end

local function transformRelaxationTime(rt, p, q, db)
   local t = {}
   t.model = rt[1]
   --
   -- V-T relaxation time models
   --
   if t.model == "Millikan-White" then
      M_p = db[p].M*1000.0 -- kg -> g
      M_q = db[q].M*1000.0 -- kg -> g
      mu = (M_p * M_q)/(M_p + M_q)
      theta_v = db[p].theta_v or db[p].vib_data.theta_v
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
      --
      -- V-V relaxation time models
      --
   elseif (t.model == "Schwartz-Slawsky-Herzfeld" or t.model == "SSH") then
      t.model = "SSH-VV"
      M_p = db[p].M
      M_q = db[q].M
      t.mu_pq = (M_p * M_q)/(M_p + M_q)
      t.mu_pp = (M_p * M_p)/(M_p + M_p)
      t.mu_qq = (M_q * M_q)/(M_q + M_q)
      t.theta_v_p = db[p].theta_v or db[p].vib_data.theta_v
      t.theta_v_q = db[q].theta_v or db[q].vib_data.theta_v
      t.sigma = 0.5*(SSHSigma(p, db) + SSHSigma(q, db))
      t.epsilon = sqrt(db[p].epsilon * db[q].epsilon)
      t.f_m_p = computeMassFactor(p, db)
      t.f_m_q = computeMassFactor(q, db)
      t.r_eq_p = db[p].r_eq
      t.r_eq_q = db[q].r_eq
   elseif t.model == "SSH-VT" then
      t.model = "SSH-VT"
      M_p = db[p].M
      M_q = db[q].M
      t.mu_pq = (M_p * M_q)/(M_p + M_q)
      t.mu_pp = (M_p * M_p)/(M_p + M_p)
      t.mu_qq = (M_q * M_q)/(M_q + M_q)
      t.theta_v_p = db[p].theta_v or db[p].vib_data.theta_v
      t.sigma = 0.5*(SSHSigma(p, db) + SSHSigma(q, db))
      t.epsilon = sqrt(db[p].epsilon * db[q].epsilon)
      t.f_m_p = computeMassFactor(p, db)
      t.r_eq_p = db[p].r_eq
   elseif t.model == "ParkN2e-" then
      t.a_low = rt.a_low
      t.b_low = rt.b_low
      t.c_low = rt.c_low
      t.a_high = rt.a_high
      t.b_high = rt.b_high
      t.c_high = rt.c_high
   elseif t.model == "Candler" then
      t.submodel = transformRelaxationTime(rt.submodel, p, q, db)
   else
      print("The relaxation time model: ", t.model, " it not known.")
      print("Bailing out!")
      os.exit(1)
   end

   return t
end

function mechanism.calculateDissociationEnergy(dissociating_species, db)
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

local function buildChemistryCouplingModel(m, participating_species, db)
   -- Coupling models for chemistry vibrational coupling need their parameters
   -- assembled into a table. We do that here, computing any things that they need.
   --
   -- Notes: Consider 
   -- @author: Nick Gibbons
   local ccm = {}
   ccm.model = m.model
   if m.model == "ImpartialDissociation" then
       ccm.Thetav = db[participating_species].theta_v or db[participating_species].vib_data.theta_v
       ccm.D = mechanism.calculateDissociationEnergy(participating_species, db)
   elseif m.model == "MarroneTreanorDissociation" then
       Runiv = 8.3145
       ccm.Thetav = m.Thetav or db[participating_species].theta_v or db[participating_species].vib_data.theta_v
      ccm.D = m.D or mechanism.calculateDissociationEnergy(participating_species, db)
      if m.U then
         ccm.U = m.U
      elseif m.TD_on_U then
         ccm.U = ccm.D / m.TD_on_U / Runiv
      else
         print("Either U or TD_on_U needs to be provided for Marrone-Treanor chemistry exchange coupling")
         print("Bailing out!")
         os.exit(1)
      end
   elseif m.model == "ModifiedMarroneTreanorDissociation" then
       Runiv = 8.3145
      ccm.Thetav = m.Thetav
      ccm.Ustar = m.Ustar
      ccm.aU = m.aU
      ccm.T_D = m.T_D or mechanism.calculateDissociationEnergy(participating_species, db) / Runiv
   elseif m.model == "HeavyParticleCollisionIonisation" then
      -- nothing extra to do
   elseif m.model == "ImpartialChem" then
      -- nothing extra to do
   elseif m.model == "ElectronImpactIonisation" then
      -- the ion is created with the ionisation energy of the species
      -- it was created from
      ps = string.gsub(participating_species, "+", "", 1)
      ccm.ionisation_energy = db[ps].ionisation_energy
   else
      print(string.format("ERROR: chemistry coupling model '%s' is not known.", m.model))
      os.exit(1)
   end
   return ccm
end

local function relaxationTimeToLuaStr(rt)
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
   elseif rt.model == "SSH-VV" then
      str = "{model='SSH-VV',\n"
      str = str .. string.format("      theta_v_p = %f,\n", rt.theta_v_p)
      str = str .. string.format("      theta_v_q = %f,\n", rt.theta_v_q)
      str = str .. string.format("      mu_pq = %.6e,\n", rt.mu_pq)
      str = str .. string.format("      mu_pp = %.6e,\n", rt.mu_pp)
      str = str .. string.format("      mu_qq = %.6e,\n", rt.mu_qq)
      str = str .. string.format("      sigma = %.6e,\n", rt.sigma)
      str = str .. string.format("      epsilon = %.6e,\n", rt.epsilon)
      str = str .. string.format("      f_m_p = %.6e,\n", rt.f_m_p)
      str = str .. string.format("      f_m_q = %.6e,\n", rt.f_m_q)
      str = str .. string.format("      r_eq_p = %.8e,\n", rt.r_eq_p)
      str = str .. string.format("      r_eq_q = %.8e,\n", rt.r_eq_q)
      str = str .. "  },"
   elseif rt.model == "SSH-VT" then
      str = "{model='SSH_VT',\n"
      str = str .. string.format("      theta_v_p = %f,\n", rt.theta_v_p)
      str = str .. string.format("      mu_pq = %.6e,\n", rt.mu_pq)
      str = str .. string.format("      mu_pp = %.6e,\n", rt.mu_pp)
      str = str .. string.format("      mu_qq = %.6e,\n", rt.mu_qq)
      str = str .. string.format("      sigma = %.6e,\n", rt.sigma)
      str = str .. string.format("      epsilon = %.6e,\n", rt.epsilon)
      str = str .. string.format("      f_m_p = %.6e,\n", rt.f_m_p)
      str = str .. string.format("      r_eq_p = %.8e,\n", rt.r_eq_p)
      str = str .. "  },"
   elseif rt.model == "Candler" then
      submodelstr = relaxationTimeToLuaStr(rt.submodel)
      str = string.format("{model='Candler', %s }", submodelstr)
   elseif rt.model == "ParkN2e-" then
      str = string.format("a_low=%.4e, b_low=%.4e, c_low=%.4e, a_high=%.4e, b_high=%.4e, c_high=%.4e", rt.a_low, rt.b_low, rt.c_low, rt.a_high, rt.b_high, rt.c_high)
   else
      print(string.format("ERROR: relaxation time model '%s' is not known.", rt.model))
      os.exit(1)
   end
   return str
end

local function chemistryCouplingTypeToLuaStr(ct)
   local str = ""
   if ct.model == "ImpartialDissociation" then
      str = string.format("{model='ImpartialDissociation', D=%f, Thetav=%f}", ct.D, ct.Thetav)
   elseif ct.model == "MarroneTreanorDissociation" then
      str = string.format("{model='MarroneTreanorDissociation', D=%f, Thetav=%f, U=%f}", ct.D, ct.Thetav, ct.U)
   elseif ct.model == "ModifiedMarroneTreanorDissociation" then
      str = string.format("{model='ModifiedMarroneTreanorDissociation', T_D=%f, Thetav=%f, Ustar=%f, aU=%f}", ct.T_D, ct.Thetav, ct.Ustar, ct.aU)
   elseif ct.model == "ImpartialChem" or ct.model == "HeavyParticleCollisionIonisation" then
      str = "{model = 'ImpartialChem'}"
   else
      print(string.format("ERROR: Chemistry coupling type '%s' is not known.", ct.model))
      os.exit(1)
   end
   return str
end

local function tableToString(o)
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

function mechanism.mechanismToLuaStr(index, m)
   local argStr
   if m.type == "V-T" or m.type == "E-V" or m.type == "V-E" then
      argStr = string.format("  p = '%s', q = '%s',\n", m.p, m.q)
      argStr = argStr .. string.format("  mode_p = %d,\n", energy_modes[m.p].vib)
      argStr = argStr .. string.format("  rate = '%s',\n", m.rate)
      argStr = argStr .. string.format("  relaxation_time = %s\n", relaxationTimeToLuaStr(m.rt))
   elseif m.type == "V-V" then
      argStr = string.format("  p = '%s', q = '%s',\n", m.p, m.q)
      argStr = argStr .. string.format("  mode_p = %d, mode_q= %d,\n", energy_modes[m.p].vib, energy_modes[m.q].vib)
      argStr = argStr .. string.format("  theta_D_p = %.3f, theta_D_q = %.3f,\n", db[m.p].vib_data.theta_D, db[m.q].vib_data.theta_D)
      argStr = argStr .. string.format("  theta_v_p = %.3f, theta_v_q = %.3f,\n", db[m.p].vib_data.theta_v, db[m.q].vib_data.theta_v)
      local R_univ = 8.31451
      argStr = argStr .. string.format("  R_p = %.3f, R_q = %.3f,\n", R_univ/db[m.p].M, R_univ/db[m.q].M)
      argStr = argStr .. string.format("  rate = '%s',\n", m.rate)
      argStr = argStr .. string.format("  relaxation_time = %s\n", relaxationTimeToLuaStr(m.rt))
   elseif m.type == "E-T" then
      argStr = string.format("  p = '%s', q = '%s',\n", m.p, m.q)
      argStr = argStr .. string.format("  rate = 'ElectronExchange',\n")
      argStr = argStr .. string.format("  mode_p = %d,\n", energy_modes[m.p].electron)
      argStr = argStr .. string.format("  exchange_cross_section = %s\n", tableToString(m.exchange_cross_section))
   elseif m.type == "C-V" or m.type == "C-E" then
      argStr = string.format("  p = '%s',\n", m.p)
      argStr = argStr .. string.format("  rate = '%s',\n", m.rate)
      argStr = argStr .. string.format("  mode_p = %d,\n", energy_modes[m.p].vib) 
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

local function parseMechString(s)
   t = lpeg.match(G, s)
   return t
end

function mechanism.validateMechanism(m, i)
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

local function expandKeywords(keyword, species, db)
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

local function expandColliders(t, species, db)
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

function mechanism.addUserMechToTable(index, m, mechanisms, species, db)
   t = parseMechString(m[1])
   ps = expandColliders(t[1], species, db)
   qs = expandColliders(t[2], species, db)

   for _,p in ipairs(ps) do
      for __,q in ipairs(qs) do
	 if not (m.type == 'V-V' and p == q) then
	    mechanisms[index] = {}
	    mechanisms[index].type = m.type
	    mechanisms[index].p = p
	    mechanisms[index].q = q
	    if m.type == 'V-T' or m.type == 'V-V' or m.type == "E-V" then 
	       mechanisms[index].rate = m.rate
	       mechanisms[index].rt = transformRelaxationTime(m.relaxation_time, p, q, db)
	    elseif m.type == 'E-T' then
	       mechanisms[index].exchange_cross_section = m.exchange_cross_section
	    end
	    index = index + 1
         end
      end
   end
   return index
end

local function seekEntriesWithLabel(table, keyname, keyvalue)
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

local function seekIons(table)
   entries = {}
   for k,v in pairs(table) do
      if string.find(k, "+") then
         entries[k]=v
      end
   end
   entries["e-"]=nil
   return entries
end

local function concatenateTables(t1, t2)
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

function mechanism.addUserChemMechToTable(index, m, mechanisms, species, db, reaction)
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

         if m.rate == "Marrone-Treanor" then
            -- for backwards compatibility, if no type was given assume C-V
            if m.type == nil then
               m.type = "C-V"
            end
            if m.type == "C-V" then
               participants = seekEntriesWithLabel(participant_species, "type", "molecule")
            else
               if m.coupling_model.model == "ElectronImpactIonisation" then
                  participants = seekIons(participant_species)
               end
               if m.coupling_model.model == "HeavyParticleCollisionIonisation" then
                   participants = seekEntriesWithLabel(participant_species, "type", "electron")
               else
                  participants = participant_species
               end
            end
         else
            print(string.format('ERROR: Unknown chemistry coupling rate: %s', m.type))
         end

         for speciesname,dbentry in pairs(participants) do
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

return mechanism
