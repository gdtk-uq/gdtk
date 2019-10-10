-- vortex-flow-spec.lua
--
-- This defines the inviscid flow field of a compressible vortex,
-- following the description in C.-W. Shu's 1996 ICASE Report 97-65
-- but using dimensional quantities.
-- Peter J. and Lachlan Whyborn 2019-10-10
--
local Rgas = 287.1    -- J/kg.K
local g    = 1.4      -- ratio of specific heats
local p0 = 100.0e3    -- Pascals
local T0 = 300.0      -- Kelvin
local rho0 = p0/(Rgas*T0) -- density, kg/m^3
local a0 = math.sqrt(g*Rgas*T0) -- sound speed, m/s
local eps = 5*a0    -- strength of vortex
local pi24 = 4*(math.pi)^2
local x_centre = 5.0
local y_centre = 5.0
L = 10.0 -- size of domain edge, in metres (global variable)

function vortex_flow(r)
   local r2 = r*r
   local tmp0 = eps^2/pi24*math.exp(1-r2)
   local v2 = tmp0*r2
   local v = math.sqrt(v2)
   local dT = (g-1)/(Rgas*g)*0.5*tmp0
   local T = T0 - dT
   local tmp1 = math.pow(rho0, g)/p0 * Rgas * T
   local rho = math.pow(tmp1, 1/(g-1) )
   local p = rho*Rgas*T
   return v, p, T, rho
end

print('Reference quantities:')
v_at1, _, _, _ = vortex_flow(1)
tau = 1/a0 
print('a0=', a0, 'eps=', eps, 'v(r=1)=', v_at1, 'tau=', tau)

-- Bulk advection velocity
velx_inf = a0; vely_inf = a0

function fillTable(tbl, x, y)
   local xbar = x - x_centre
   local ybar = y - y_centre
   local r = math.sqrt(xbar^2 + ybar^2)
   local speed, press, TKelvin, density = vortex_flow(r)
   tbl.p = press
   if r < 1.0e-9 then
      tbl.velx = velx_inf
      tbl.vely = vely_inf
   else
      tbl.velx = -ybar/r * speed + velx_inf
      tbl.vely = xbar/r * speed + vely_inf
   end
   tbl.velz = 0.0
   tbl.T = TKelvin
   -- omit mass fractions; we have a single species gas model
   return tbl
end


