local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function fillTable(tab, x, y, t)
   

tab.p = 50000.0*sin(3.1415926535897932*y) + 20000.0*cos(6.2831853071795865*x) + 100000.0



tab.T = 0.00348432055749129*(50000.0*sin(3.1415926535897932*y) +20000.0*cos(6.2831853071795865*x) + 100000.0)/(0.15*sin(3.1415926535897932*x) - 0.1*cos(1.5707963267948966*y) +1.0)



tab.velx = 50.0*sin(4.7123889803846899*x) - 30.0*cos(1.8849555921538759*y) + 800.0



tab.vely = 40.0*sin(2.0943951023931954*y) - 75.0*cos(1.5707963267948966*x) + 800.0



   tab.velz = 0.0
   return tab
end

function ghostCells(args)
   ghost0 = {}
   fillTable(ghost0, args.gc0x, args.gc0y,t)
   return ghost0
end

function interface(args)
   x = args.x; y = args.y; t = args.t
   iface = {}
   fillTable(iface, x, y, t);
   return iface
end


