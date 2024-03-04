local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function solidInterface(args)
   x = args.x; y = args.y
   

T_s = 2.5*sin(3.1415926535897932*y)*cos(2.3561944901923449*x) -10.0*cos(4.7123889803846899*x) + 350



   return T_s
end

