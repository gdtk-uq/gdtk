-- udf-supersonic-in.lua
-- Lua script for the user-defined functions called
-- by the UserDefinedBC in the mallison.lua script.
--
-- This script allows us to specify a supersonic inflow
-- condition which varies across the inflow boundary,
-- and which is defined by a function.
--
-- Wilson Chan 05-Nov-2008


function ghostCells(args)
   -- Some brief calculations to obtain inflow conditions at different
   -- axial locations - Start at cylinder leading edge at Mach 8.8 with 
   -- an increasing axial mach no. gradient of 0.24/m (as specified in 
   -- papers by Mallinson et al. and Boyce & Hillier).
   mach_inf = 8.8 + 0.24 * args.x
   -- Isentropic calculation, T_inf = T_0 / (1 + (g-1)/2*M^2)
   T_inf = 1150.0 / (1.0 + (1.4-1.0)/2.0 * mach_inf^2)
   -- Isentropic calculation, p_inf = p_0 / (T_0/T)^(g/(g-1))
   p_inf = 60.0e6 / (1150.0/T_inf)^3.5
   -- Compute velocity from temperature and Mach number
   vel_nominal = mach_inf * (1.4 * 297 * T_inf)^0.5
   --
   -- Turbulence quantities for free stream 
   tke_inf = 1.0e-12
   omega_inf = 1.0   
   --
   -- Input inflow conditions for each cell along boundary
   ghost = {}
   ghost.p = p_inf -- pressure, Pa
   ghost.T = T_inf
   ghost.velx = vel_nominal * math.cos(math.rad(5.0 * args.y))
   ghost.vely = vel_nominal * math.sin(math.rad(5.0 * args.y))
   ghost.velz = 0.0
   ghost.tke = tke_inf
   ghost.omega = omega_inf
   ghost.massf = {N2=1.0} -- mass fractions to be provided as a table
   return ghost, ghost
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- Some brief calculations to obtain inflow conditions at different
   -- axial locations - Start at cylinder leading edge at Mach 8.8 with 
   -- an increasing axial mach no. gradient of 0.24/m (as specified in 
   -- papers by Mallinson et al. and Boyce & Hillier).
   mach_inf = 8.8 + 0.24 * args.x
   -- Isentropic calculation, T_inf = T_0 / (1 + (g-1)/2*M^2)
   T_inf = 1150.0 / (1.0 + (1.4-1.0)/2.0 * mach_inf^2)
   vel_nominal = mach_inf * (1.4 * 297 * T_inf)^0.5
   --
   -- Turbulence quantities for free stream 
   tke_inf = 1.0e-12
   omega_inf = 1.0
   --
   -- Input inflow conditions for each interface along boundary
   face = {}
   face.velx = vel_nominal * math.cos(math.rad(5.0 * args.y))
   face.vely = vel_nominal * math.sin(math.rad(5.0 * args.y))
   face.velz = 0.0
   face.tke = tke_inf
   face.omega = omega_inf
   face.T = T_inf
   face.massf = {N2=1.0} -- mass fractions to be provided as a table
   return face
end
