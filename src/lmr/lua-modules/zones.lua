-- Module for constructing special zone objects, required by prep.lua.
--
-- Authors: PJ and RJG
--

-- Classes for construction of zones.

local ReactionZone = {
   p0 = nil,
   p1 = nil
}

function ReactionZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to ReactionZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(reactionZones)
   reactionZones[#(reactionZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

local IgnitionZone = {
   p0 = nil,
   p1 = nil,
   T = nil -- degrees K
}

function IgnitionZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1", "T"})
   if not flag then
      error("Invalid name for item supplied to IgnitionZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(ignitionZones)
   ignitionZones[#(ignitionZones)+1] = o
   -- Must have corners and temperature
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   if not o.T then
      error("You need to supply ignition temperature T", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

local TurbulentZone = {
   p0 = nil,
   p1 = nil,
}

function TurbulentZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to TurbulentZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(turbulentZones)
   turbulentZones[#(turbulentZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

local SuppressReconstructionZone = {
   p0 = nil,
   p1 = nil,
}

function SuppressReconstructionZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to SuppressReconstructionZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(suppressReconstructionZones)
   suppressReconstructionZones[#(suppressReconstructionZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

local SuppressViscousStressesZone = {
   p0 = nil,
   p1 = nil,
}

function SuppressViscousStressesZone:new(o)
   o = o or {}
   local flag = checkAllowedNames(o, {"p0", "p1"})
   if not flag then
      error("Invalid name for item supplied to SuppressViscousStressesZone constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new zone, for later construction of the config file.
   -- Note that we want zone id to start at zero for the D code.
   o.id = #(suppressViscousStressesZones)
   suppressViscousStressesZones[#(suppressViscousStressesZones)+1] = o
   -- Must have corners
   if not o.p0 then
      error("You need to supply lower-left corner p0", 2)
   end
   if not o.p1 then
      error("You need to supply upper-right corner p1", 2)
   end
   o.p0.z = o.p0.z or 0.0
   o.p1.z = o.p1.z or 0.0
   return o
end

return {
   ReactionZone = ReactionZone,
   IgnitionZone = IgnitionZone,
   TurbulentZone = TurbulentZone,
   SuppressReconstructionZone = SuppressReconstructionZone,
   SuppressViscousStressesZone = SuppressViscousStressesZone
}
