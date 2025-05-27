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


fmass = -1000.0*pi*sin(6283.1853071795865*t)



fxmom = -180000.0*pi*(0.5*cos(6283.1853071795865*t) + 1.0)^2*sin(6283.1853071795865*t)



fymom = -30000.0*pi*(0.5*sin(6283.1853071795865*t) + 1.0)^2*sin(6283.1853071795865*t) + 60000.0*pi*(0.5*sin(6283.1853071795865*t) + 1.0)*(0.5*cos(6283.1853071795865*t) + 1.0)*cos(6283.1853071795865*t)



fe = (0.5*cos(6283.1853071795865*t) + 1.0)*(1800000.0*pi*(0.5*sin(6283.1853071795865*t) + 1.0)^3*cos(6283.1853071795865*t) + 2500.0*pi*(50000.0*sin(6283.1853071795865*t) +100000.0)*sin(6283.1853071795865*t)/(0.5*cos(6283.1853071795865*t) + 1.0)^2 - 7200000.0*pi*(0.5*cos(6283.1853071795865*t) + 1.0)^3*sin(6283.1853071795865*t) +250000000.0*pi*cos(6283.1853071795865*t)/(0.5*cos(6283.1853071795865*t) + 1.0)) - 1000.0*pi*(450.0*(0.5*sin(6283.1853071795865*t) + 1.0)^4 + 2.5*(50000.0*sin(6283.1853071795865*t) + 100000.0)/(0.5*cos(6283.1853071795865*t) + 1.0) + 1800.0*(0.5*cos(6283.1853071795865*t) + 1.0)^4)*sin(6283.1853071795865*t)



   src.mass = fmass
   src.momentum_x = fxmom
   src.momentum_y = fymom
   src.momentum_z = 0.0
   src.total_energy = fe
   return src
end
