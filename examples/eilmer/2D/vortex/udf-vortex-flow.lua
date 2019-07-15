-- udf-vortex-flow.lua
-- Lua script for the inviscid vortex flow.
--
-- This particular example defines the inviscid flow field
-- of a compressible vortex.

Rgas  = 287      -- J/kg.K
g     = 1.4      -- ratio of specific heats
-- radial limits of flow domain
r_i   = 1.0      -- metres
r_o   = 1.384
-- Set flow properties ar the inner radius.
p_i   = 100.0e3                  -- Pa
M_i   = 2.25
rho_i = 1.0                      -- kg/m**3
T_i   = p_i / (Rgas * rho_i)     -- K
a_i   = math.sqrt(g * Rgas * T_i)     -- m/s
u_i   = M_i * a_i                -- m/s

if false then
   print("Set up inviscid vortex")
   print("    p_i=", p_i, "M_i=", M_i, "rho_i=", rho_i, 
	 "T_i=", T_i, "a_i=", a_i, "u_i=", u_i)
end

function vortex_flow(r)
   local u   = u_i * r_i / r
   local t1  = r_i / r
   local t2  = 1.0 + 0.5 * (g - 1.0) * M_i * M_i * (1.0 - t1 * t1)
   local rho = rho_i * math.pow( t2, 1.0/(g - 1.0) )
   local p = p_i * math.pow( rho/rho_i, g )
   local T = p / (rho * Rgas)
   return u, p, T
end

function fillTable(t, x, y)
   local r = math.sqrt(x*x + y*y)
   local theta = math.atan2(y, x)
   local speed, press, TKelvin = vortex_flow(r)
   t.p = press
   t.velx = math.sin(theta) * speed
   t.vely = -math.cos(theta) * speed
   t.velz = 0.0
   t.T = TKelvin
   -- omit mass fractions; we have a single species gas model
   return t
end

function twoGhostCells(args)
   ghost0 = fillTable({}, args.gc0x, args.gc0y)
   ghost1 = fillTable({}, args.gc1x, args.gc1y)
   return ghost0, ghost1
end

-- The following functions are expected by the boundary-condition
-- interpreter.  It will call them up at run time, at the appropriate
-- stage of the time-step calculation.

function ghostCells(args)
   return twoGhostCells(args)
end

function interface(args)
   return fillTable({}, args.x, args.y)
end

