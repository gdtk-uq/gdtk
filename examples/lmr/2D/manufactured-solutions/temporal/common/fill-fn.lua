local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi
t = 0.0
function gasFillFn(x, y, z)
   

p = 50000.0*sin(6283.1853071795865*t) + 100000.0



T  = 0.00348432055749129*(50000.0*sin(6283.1853071795865*t) +100000.0)/(0.5*cos(6283.1853071795865*t) + 1.0)



velx = 60.0*(0.5*cos(6283.1853071795865*t) + 1.0)^2



vely = 30.0*(0.5*sin(6283.1853071795865*t) + 1.0)^2



   return FlowState:new{p=p, T=T, velx=velx, vely=vely}
end
