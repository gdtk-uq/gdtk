-- udf-source-template.lua
-- Lua template for the source terms of a Manufactured Solution.
--
-- PJ, 29-May-2011
-- RJG, 06-Jun-2014
--   Declared maths functions as local

local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

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
   return src
end
