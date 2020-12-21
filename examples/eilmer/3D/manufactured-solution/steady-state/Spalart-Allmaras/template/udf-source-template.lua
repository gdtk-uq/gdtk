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
local min = math.min
local max = math.max
local sqrt = math.sqrt
local abs = math.abs
local tanh = math.tanh

function sourceTerms(t, cell)
   src = {}
   x = cell.x
   y = cell.y
   z = cell.z
   $expressions
   src.mass = fmass
   src.momentum_x = fxmom
   src.momentum_y = fymom
   src.momentum_z = fzmom
   src.total_energy = fe
   src.nuhat = fnuhat
   return src
end
