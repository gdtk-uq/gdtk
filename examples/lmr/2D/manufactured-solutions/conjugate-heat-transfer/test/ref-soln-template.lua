local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function refSoln(t, x, y, z)
   tab = {}
   <insert-expressions-here>
   return tab
end

function refSolidSoln(t, x, y, z)
   tab = refSoln(t, x, y, z)
   -- Solid domain expects "T"
   tab.T = tab.T_s
   return tab
end

