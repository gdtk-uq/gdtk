local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi
local max = math.max
local sqrt = math.sqrt
local tanh = math.tanh

function fillTable(tab, x, y, z, t)
   $expressions
   return tab
end

function ghostCells(args)
   x = args.x; y = args.y; z = args.z; t = args.t
   i = args.i; j = args.j; k = args.k
   ghost = {}
   fillTable(ghost, args.gc0x, args.gc0y, args.gc0z, t)
   return ghost
end

function interface(args)
   x = args.x; y = args.y; z = args.z; t = args.t
   iface = {}
   fillTable(iface, x, y, z, t);
   return iface
end


