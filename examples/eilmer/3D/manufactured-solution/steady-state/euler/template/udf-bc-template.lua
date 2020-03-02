local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function fillTable(tab, x, y, z, t)
   tab.T = {}
   <insert-expressions-here>
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


