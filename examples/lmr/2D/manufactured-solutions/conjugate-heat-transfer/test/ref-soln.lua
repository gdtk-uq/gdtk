local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function refSoln(t, x, y, z)
   tab = {}
   

tab.rho = 0.1*sin(2.3561944901923449*x) + 0.15*cos(3.1415926535897932*y) + 0.08*cos(3.9269908169872415*x*y) +1.0



tab.p = (25.0*sin(3.1415926535897932*y)*cos(2.3561944901923449*x)- 10.0*cos(4.7123889803846899*x) + 350)*(28.7*sin(2.3561944901923449*x) + 43.05*cos(3.1415926535897932*y) +22.96*cos(3.9269908169872415*x*y) + 287.0)



tab.T = 25.0*sin(3.1415926535897932*y)*cos(2.3561944901923449*x) -10.0*cos(4.7123889803846899*x) + 350



tab.T_s = 2.5*sin(3.1415926535897932*y)*cos(2.3561944901923449*x)- 10.0*cos(4.7123889803846899*x) + 350



tab['vel.x'] = 0.1*cos(5.235987755982989*x) + 1.0*cos(3.1415926535897932*y) - 0.1*cos(5.235987755982989*x*y) +1.0



tab['vel.y'] = -0.9*sin(1.5707963267948966*y) - 0.02*cos(4.7123889803846899*x) + 0.02*cos(4.7123889803846899*x*y) +0.9



   return tab
end

function refSolidSoln(t, x, y, z)
   tab = refSoln(t, x, y, z)
   -- Solid domain expects "T"
   tab.T = tab.T_s
   return tab
end

