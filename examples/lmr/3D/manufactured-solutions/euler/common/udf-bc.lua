local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function fillTable(tab, x, y, z, t)
   tab.T = {}
   

tab.p = 50000.0*sin(3.1415926535897932*y) - 35000.0*sin(1.0471975511965977*z) + 20000.0*cos(6.2831853071795865*x) +100000.0



tab.T = 0.00348432055749129*(50000.0*sin(3.1415926535897932*y) -35000.0*sin(1.0471975511965977*z) + 20000.0*cos(6.2831853071795865*x) + 100000.0)/(0.15*sin(3.1415926535897932*x) - 0.1*cos(1.5707963267948966*y) -0.12*cos(4.7123889803846899*z) + 1.0)



tab.velx = 50.0*sin(4.7123889803846899*x) - 30.0*cos(1.8849555921538759*y) - 18.0*cos(1.5707963267948966*z) +800.0



tab.vely = -75.0*sin(1.5707963267948966*x) + 40.0*cos(2.0943951023931954*y) - 30.0*cos(3.9269908169872415*z) +800.0



tab.velz = 15.0*sin(1.0471975511965977*x) - 25.0*cos(4.7123889803846899*y) + 35.0*cos(3.1415926535897932*z) +800.0



   return tab
end

function ghostCells(args)
   x = args.x; y = args.y; z = args.z; t = args.t
   i = args.i; j = args.j; k = args.k
   ghost0 = {}
   fillTable(ghost0, args.gc0x, args.gc0y, args.gc0z, t)
   --ghost1 = {}
   --fillTable(ghost1, args.gc1x, args.gc1y, args.gc1z, t)
   return ghost0--, ghost1
end

function interface(args)
   x = args.x; y = args.y; z = args.z; t = args.t
   iface = {}
   fillTable(iface, x, y, z, t);
   return iface
end


