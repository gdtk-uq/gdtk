-- udf-source-template.lua
-- Lua template for the source terms of a Manufactured Solution.
--
-- PJ, 29-May-2011
-- RJG, 06-Jun-2014
--   Declared maths functions as local

local function Heaviside(x)
   if x < 0.0 then 
      return 0.0
   else
      return 1.0
   end
end

local function Min1(x, y)
   if y > x then
       return 0.125
   else
       return 0.0
   end
end

local function Min2(x, y)
   if y > x then
       return x
   else
       return y
   end
end

local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi
local max = math.max
local sqrt = math.sqrt
local min = math.min

function sourceTerms(t, cell)
   src = {}
   x = cell.x
   y = cell.y
<insert-expressions-here>
   src.mass = fmass
   src.momentum_x = fxmom
   src.momentum_y = fymom
   src.momentum_z = 0.0
   src.total_energy = fe
   src.tke = ftke
   src.omega = fomega
   return src
end
