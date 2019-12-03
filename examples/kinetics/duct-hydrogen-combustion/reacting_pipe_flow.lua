-- reacting_pipe_flow.lua
--
-- Combustion in a supersonic stream of constant cross-section.
-- This is set up to approximate the Bittker-Scullin case 3,
-- as used by Fabian Zander in the hydrogen-combustion test case. 
--
-- PJ 2011-06-21 first version written for Eilmer3 test case
--    2016-03-19 adapted from reacting-pipe-flow.py for Eilmer3
--    2018-04-21 updated to use current gas model and reaction calls
--
-- Run with the commands:
-- $ prep-gas combusting-species.inp h2-o2-n2-9sp.lua
-- $ prep-chem h2-o2-n2-9sp.lua Bittker-Scullin.lua h2-o2-n2-9sp-18r.lua
-- $ gas-calc reacting_pipe_flow.lua

-------------------------------------------------------------
-- Things that we'll make use of shortly.

sample_header = "# x(m) rho(kg/m**3) p(Pa) T(degK) e(J/kg) v(m/s) "..
   "massf_OH massf_H2O dt_suggest(s)"

function sample_data(x, v, gas, dt_suggest)
   return string.format("%g %g %g %g %g %g %g %g %g ",
			x, gas.rho, gas.p, gas.T, gas.u,
			v, gas.massf["OH"], gas.massf["H2O"], dt_suggest)
end

function eos_derivatives(gas0, gmodel, tol)
   -- Finite difference evaluation, assuming that gas0 is valid state already.
   if not tol then tol = 0.0001 end
   local gas1 = GasState:new{gmodel}
   copyValues(gas0, gas1)
   local p0 = gas0.p
   local rho0 = gas0.rho
   local u0 = gas0.u
   --
   local drho = rho0 * tol
   gas1.rho = rho0 + drho
   gmodel:updateThermoFromRHOU(gas1)
   local dpdrho = (gas1.p - p0)/drho
   --
   gas1.rho = rho0
   local du = u0 * tol
   gas1.u = u0 + du
   gmodel:updateThermoFromRHOU(gas1)
   local dpdu = (gas1.p - p0)/du
   --
   return dpdrho, dpdu
end

------------------------------------------------------------------
-- Start the main script...
debug = false
print("# Reacting pipe flow -- Bittker-Scullin test case 3.")
print(sample_header)

-- Gas model setup
gmodel = GasModel:new{"h2-o2-n2-9sp.lua"}
nsp = gmodel:nSpecies()
nmodes = gmodel:nModes()
if debug then
   print("nsp=", nsp, " nmodes=", nmodes, " gmodel=", gmodel)
end

print("# Gas properties at the start of the pipe.")
local gas0 = GasState:new{gmodel}
gas0.p = 96.87e3 -- Pa
local x = 0.0 -- m  (inlet of pipe)
local v = 4551.73 -- m/s
gas0.T = 1559.0 -- degree K
molefInit = {O2=0.1480, N2=0.5562, H2=0.2958}
gas0.massf = gmodel:molef2massf(molefInit)
gmodel:updateThermoFromPT(gas0)
local dt_suggest = 1.0e-8  -- suggested starting time-step for chemistry updater
print(sample_data(x, v, gas0, dt_suggest))

print("# Start reactions...")
local chemUpdate = ChemistryUpdate:new{filename="h2-o2-n2-9sp-18r.lua",
                                       gasmodel=gmodel}
local t = 0 -- time is in seconds
local t_final = 22.0e-6
local t_inc = 0.05e-6
local nsteps = math.floor(t_final / t_inc)
for j=1, nsteps do
   -- At the start of the step...
   local rho = gas0.rho
   local T = gas0.T
   local p = gas0.p
   local u = gas0.u
   --
   -- Do the chemical increment.
   local gas1 = GasState:new{gmodel} -- make the new one as a clone
   copyValues(gas0, gas1)
   dt_suggest = chemUpdate:updateState(gas1, t_inc, dt_suggest, gmodel)
   gmodel:updateThermoFromRHOU(gas1)
   --
   local du_chem = gas1.u - u
   local dp_chem = gas1.p - p
   if debug then print("# du_chem=", du_chem, "dp_chem=", dp_chem) end
   --
   -- Do the gas-dynamic accommodation after the chemical change.
   local Etot = u + 0.5*v*v
   local dfdr, dfdu = eos_derivatives(gas1, gmodel)
   if debug then print("# dfdr=", dfdr, "dfdu=", dfdu) end
   --[=[ Linear solve to get the accommodation increments.
      [v,      rho,        0.0, 0.0  ]   [drho  ]   [0.0           ]
      [0.0,    rho*v,      1.0, 0.0  ] * [dv    ] = [-dp_chem      ]
      [v*Etot, rho*Etot+p, 0.0, rho*v]   [dp_gda]   [-rho*v*du_chem]
      [dfdr,   0.0,       -1.0, dfdu ]   [du_gda]   [0.0           ]
   ]=]
   -- Compute the accommodation increments using expressions from Maxima.
   local denom = rho*rho*v*v - dfdr*rho*rho - dfdu*p
   local drho = (dp_chem - du_chem*dfdu)*rho*rho / denom
   local dv = -(dp_chem - du_chem*dfdu)*rho*v / denom
   local dp_gda = -(du_chem*dfdu*rho*rho*v*v - dfdr*dp_chem*rho*rho 
                       - dfdu*dp_chem*p) / denom
   local du_gda = -(du_chem*rho*rho*v*v - du_chem*dfdr*rho*rho - dp_chem*p) / denom
   if debug then 
      print("# drho=", drho, "dv=", dv, "dp_gda=", dp_gda, "du_gda=", du_gda)
      print("# residuals=", v*drho + rho*dv, rho*v*dv + dp_gda + dp_chem,
            v*Etot*drho + (rho*Etot+p)*dv + rho*v*du_gda + rho*v*du_chem,
            dfdr*drho - dp_gda + dfdu*du_gda)
   end
   -- Add the accommodation increments.
   gas1.rho = gas0.rho + drho
   v1 = v + dv
   p1_check = gas1.p + dp_gda
   gas1.u = gas1.u + du_gda
   gmodel:updateThermoFromRHOU(gas1)
   if debug then
      print("# At new point for step ", j, ": gas1.p=", gas1.p, "p1_check=", p1_check, 
            "rel_error=", math.abs(gas1.p-p1_check)/p1_check)
   end
   -- Have now finished the chemical and gas-dynamic update.
   t = t + t_inc
   x = x + 0.5*(v + v1) * t_inc
   print(sample_data(x, v1, gas1, dt_suggest))
   -- House-keeping for the next step.
   v = v1
   copyValues(gas1, gas0) -- gas0 will be used in the next iteration
end
print("# Done.")
