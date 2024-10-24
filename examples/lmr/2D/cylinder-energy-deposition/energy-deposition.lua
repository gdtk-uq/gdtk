-- energy-deposition.lua
dofile('parameters.lua')

Q = 4.0e3 -- W for a 1m depth flow
xC = -0.9*h
xL = xC - (0.5*h/15)
xR = xC + (0.5*h/15)
yE = 0.5*h/56
q_dot = Q/(h/15)/(h/56) -- rate of heat added per volume

function sourceTerms(t, cell)
   local src = {}
   local x = cell.x
   local y = cell.y
   src.mass = 0.0
   src.momentum_x = 0.0
   src.momentum_y = 0.0
   src.total_energy = 0.0
   if (y < yE) and (x >= xL) and (x <= xR) then
      src.total_energy = q_dot
   end
   return src
end
