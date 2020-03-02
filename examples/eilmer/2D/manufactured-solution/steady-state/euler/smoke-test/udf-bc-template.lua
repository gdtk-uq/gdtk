local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi

function fillTable(tab, x, y, t)
   tab.T = {}
   <insert-expressions-here>
   tab.velz = 0.0
   return tab
end

function ghostCells(args)
   x = args.x; y = args.y; t = args.t
   i = args.i; j = args.j; k = args.k
   ghost0 = {}
   fillTable(ghost0, args.gc0x, args.gc0y,t)
   --ghost1 = {}
   --fillTable(ghost1, args.gc1x, args.gc1y,t)
   return ghost0--, ghost1
end

function interface(args)
   x = args.x; y = args.y; t = args.t
   iface = {}
   fillTable(iface, x, y, t);
   return iface
end


