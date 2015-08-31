local pi = math.pi
local cos = math.cos
local sin = math.sin

local L = 1.0
local R = 287.0
local gam = 1.4
local k_g = 10000.0
local k_s = 10*k_g

rho0=1.0; rhox=0.1; rhoy=0.15; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
u0 = 1.0; ux = 0.1; uy = u0; uxy = -ux; aux = 5.0/3; auy = -1.0; auxy = aux;
v0 = 0.9; vx = -0.02; vy = -v0; vxy = -vx; avx = 1.5; avy = 0.5; avxy = avx;	
T0 = 350; Tx = -10.0; Ty = 25.0; aTx = 1.5; aTy = 1.0; Ti = 350.0; aTx2 = 0.75;

function rho(x, y)
   return rho0 + rhox*sin(arhox*pi*x/L) + rhoy*cos(arhoy*pi*y/L) + rhoxy*cos(arhoxy*pi*x*y/(L*L));
end

function u(x, y)
   return u0 + ux*cos(aux*pi*x/L) + uy*cos(auy*pi*y/L) + uxy*cos(auxy*pi*x*y/(L*L))
end

function v(x, y)
   return  v0 + vx*cos(avx*pi*x/L) + vy*sin(avy*pi*y/L) + vxy*cos(avxy*pi*x*y/(L*L))
end

function T(x, y)
   return T0 + Tx*cos(aTx*pi*x/L) + Ty*cos(aTx2*pi*x/L)*sin(aTy*pi*y/L)
end

function Ts(x, y)
   return T0 + Tx*cos(aTx*pi*x/L) +  Ty*(k_g/k_s)*cos(aTx2*pi*x/L)*sin(aTy*pi*y/L)
end

function refSoln(x, y, z)
   t = {}
   t["T[0]"] = T(x, y)
   return t
end

function refSolidSoln(x, y, z)
   t = {}
   t.T = Ts(x, y)
   return t
end

function fillTable(t, x, y)
   t_rho = rho(x, y)
   t.T = {}
   t.T[1] = T(x,y)
   t.p = t_rho*R*t.T[1]
   t.velx = u(x, y)
   t.vely = v(x, y)
   t.velz = 0.0
   t.massf = {1.0}     -- mass fractions to be provided as a table
   return t
end

function ghostCells_north(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains {t, x, y, z, csX, csY, csZ, i, j, k}
   -- Set constant conditions across the whole boundary.
   x = args.x; y = args.y
   i = args.i; j = args.j; k = args.k
   ghost1 = {}
   ghost2 = {}
   cell = sampleFlow(blkId, i, j+1, k)
   ghost1 = fillTable(ghost1, cell.x, cell.y)
   cell = sampleFlow(blkId, i, j+2, k)
   ghost2 = fillTable(ghost2, cell.x, cell.y)
   return ghost1, ghost2
end

function ghostCells_east(args)
   x = args.x; y = args.y
   i = args.i; j = args.j; k = args.k
   ghost1 = {}
   ghost2 = {}
   cell = sampleFlow(blkId, i+1, j, k)
   ghost1 = fillTable(ghost1, cell.x, cell.y)
   cell = sampleFlow(blkId, i+2, j, k)
   ghost2 = fillTable(ghost2, cell.x, cell.y)
   return ghost1, ghost2
end

function ghostCells_south(args)
   x = args.x; y = args.y
   i = args.i; j = args.j; k = args.k
   ghost1 = {}
   ghost2 = {}
   cell = sampleFlow(blkId, i, j-1, k)
   ghost1 = fillTable(ghost1, cell.x, cell.y)
   cell = sampleFlow(blkId, i, j-2, k)
   ghost2 = fillTable(ghost2, cell.x, cell.y)
   return ghost1, ghost2
end

function ghostCells_west(args)
   x = args.x; y = args.y
   i = args.i; j = args.j; k = args.k
   ghost1 = {}
   ghost2 = {}
   cell = sampleFlow(blkId, i-1, j, k)
   ghost1 = fillTable(ghost1, cell.x, cell.y)
   cell = sampleFlow(blkId, i-2, j, k)
   ghost2 = fillTable(ghost2, cell.x, cell.y)
   return ghost1, ghost2
end

function interface(args)
   x = args.x; y = args.y;
   iface = {}
   fillTable(iface, x, y);
   return iface
end

interface_north = interface
interface_east = interface
interface_south = interface
interface_west = interface

function solidInterface(args)
   return Ts(args.x, args.y)
end

solidInterface_north = solidInterface
solidInterface_east = solidInterface
solidInterface_south = solidInterface
solidInterface_west = solidInterface
