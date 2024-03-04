local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function solidSourceTerms(t, cell)
   src = {}
   x = cell.x
   y = cell.y


fe_s = 390625.0*pi^2*sin(3.1415926535897932*y)*cos(2.3561944901923449*x) - 2250000.0*pi^2*cos(4.7123889803846899*x)



   return fe_s
end
