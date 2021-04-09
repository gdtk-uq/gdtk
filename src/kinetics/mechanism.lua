-- Author: Rowan J. Gollan
-- Date: 2021-04-08
--

module(..., package.seeall)

local pow = math.pow
local sqrt = math.sqrt

require 'lex_elems'

function transformRelaxationTime(rt, p, q, db)
   t = {}
   t.model = rt[1]
   if t.model == "Millikan-White" then
      M_p = db[p].M*1000.0 -- kg -> g
      M_q = db[q].M*1000.0 -- kg -> g
      mu = (M_p * M_q)/(M_p + M_q)
      theta_v = db[p].theta_v
      t.a = 1.16e-3*sqrt(mu)*pow(theta_v, 4/3)
      t.b = 0.015*pow(mu, 1/4)
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
   else
      print(string.format("ERROR: relaxation time model '%s' is not known.", rt.model))
      os.exit(1)
   end
   return str
end

function mechanismToLuaStr(m, p, q)
   local typeStr
   if m.type == "V-T" then
      typeStr = "VT"
   else
      print("ERROR: type is not known: ", type)
      os.exit(1)
   end
   local key = p .. ":" .. q .. "|" .. typeStr
   local mstr = string.format("mechanism['%s'] = {\n", key)
   mstr = mstr .. string.format("  rate = '%s',\n", m.rate)
   mstr = mstr .. string.format("  relaxation_time = %s\n", relaxationTimeToLuaStr(m.rt))
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

function expandPColliders(t, molecules)
   ps = {}
   unique_ps = {}
   for i,p in ipairs(t) do
      if p == "*molcs" then
         for _,m in ipairs(molecules) do
            if unique_ps[m] == nil then
               unique_ps[m] = true
               ps[#ps+1] = m
            end
         end
      else
         if unique_ps[p] == nil then
            unique_ps[p] = true
            ps[#ps+1] = p
         end
      end
   end
   return ps
end

function expandQColliders(t, species)
   qs = {}
   unique_qs = {}
   for i,q in ipairs(t) do
      if q == "*all" then
         for _,s in ipairs(species) do
            if unique_qs[s] == nil then
               unique_ps[s] = true
               qs[#qs+1] = s
            end
         end
      else
         if unique_qs[q] == nil then
            unique_qs[q] = true
            qs[#qs+1] = q
         end
      end
   end
   return qs
end

function addUserMechToTable(m, mechanisms, species, molecules, db)
   t = parseMechString(m[1])
   ps = expandPColliders(t[1], molecules)
   qs = expandQColliders(t[2], species)

   for _,p in ipairs(ps) do
      if not mechanisms[p] then
         mechanisms[p] = {}
      end
      
      for __,q in ipairs(qs) do
         mechanisms[p][q] = {}
         mechanisms[p][q].type = m.type
         mechanisms[p][q].rate = m.rate
         mechanisms[p][q].rt = transformRelaxationTime(m.relaxation_time, p, q, db)
      end
   end
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
