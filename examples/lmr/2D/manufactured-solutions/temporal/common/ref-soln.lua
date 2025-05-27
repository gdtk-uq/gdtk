local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function refSoln(t, x, y, z)
   tab = {}
   

tab.rho = 0.5*cos(6283.1853071795865*t) + 1.0



tab.p = 50000.0*sin(6283.1853071795865*t) + 100000.0



tab.T = 0.00348432055749129*(50000.0*sin(6283.1853071795865*t) +100000.0)/(0.5*cos(6283.1853071795865*t) + 1.0)



tab['vel.x'] = 60.0*(0.5*cos(6283.1853071795865*t) + 1.0)^2



tab['vel.y'] = 30.0*(0.5*sin(6283.1853071795865*t) + 1.0)^2



   return tab
end

refSolidSoln = refSoln

