-- udf-source-template.lua
-- Lua template for the source terms of a Manufactured Solution.
--
-- PJ, 29-May-2011
-- RJG, 06-Jun-2014, Declared maths functions as local
-- PJ, 2019-10-09, sample at multiple interior points

local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

local multipoint_sample = true

function sourceTermsAtPoint(x, y)
<insert-expressions-here>
      return fmass, fxmom, fymom, fe
end

function sourceTerms(t, cell)
   src = {}
   if multipoint_sample then
      -- Sample at two Gauss points for case 6 1D high-order test
      local dx = cell.iLength / (2*math.sqrt(3))
      fmassA, fxmomA, fymomA, feA = sourceTermsAtPoint(cell.x-dx, cell.y)
      fmassB, fxmomB, fymomB, feB = sourceTermsAtPoint(cell.x+dx, cell.y)
      src.mass = 0.5*(fmassA + fmassB)
      src.momentum_x = 0.5*(fxmomA + fxmomB)
      src.momentum_y = 0.5*(fymomA + fymomB)
      src.momentum_z = 0.0
      src.total_energy = 0.5*(feA + feB)
   else
      -- Single sample at midpoint sample.
      fmass, fxmom, fymom, fe = sourceTermsAtPoint(cell.x, cell.y)
      src.mass = fmass
      src.momentum_x = fxmom
      src.momentum_y = fymom
      src.momentum_z = 0.0
      src.total_energy = fe
   end
   return src
end
