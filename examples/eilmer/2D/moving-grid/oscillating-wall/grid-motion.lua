-- Vertex velocities for oscillating plate.
-- Peter J. 2022-02-26
--
local niv = 31 -- Needs to match the number in the input script.
local vely = {}; for i=1,niv do vely[i] = 0.0 end
local Vmax = 20.0
local fHz = 100.0
local w = 2.0*math.pi*fHz

function assignVtxVelocities(t, dt)
   for i=1,niv do
      vely[i] = Vmax * math.sin((i-1)/(niv-1)*math.pi) * math.sin(w*t)
   end
   setVtxVelocitiesForBlockBoundary(1, "south", {}, vely, {})
end
