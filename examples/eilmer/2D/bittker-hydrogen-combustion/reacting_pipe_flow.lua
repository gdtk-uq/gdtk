-- reacting_pipe_flow.lua
--
-- Combustion in a supersonic stream of constant cross-section.
-- This is set up to approximate the Bittker-Scullin case 3
-- as used by Fabs in the hydrogen-combustion test case. 
--
-- PJ 2016-03-19 adapted from reacting-pipe-flow.py for Eilmer3
--
-- Invoke with the command line:
-- $ e4shared --custom-post --script-file=reacting_pipe_flow.lua

-------------------------------------------------------------
-- Things that we'll make use of shortly.

sample_header = "# x(m) rho(kg/m**3) p(Pa) T(degK) e(J/kg) v(m/s) "..
   "massf_OH massf_H2O dt_suggest(s)"

function sample_data(x, v, gas, dt_suggest)
   return string.format("%g %g %g %g %g %g %g %g %g ",
			x, gas.rho, gas.p, gas.T[1], gas.e[1],
			v, gas.massf["OH"], gas.massf["H2O"], dt_suggest)
end

function eos_derivatives(gas0, gmodel, tol)
   -- Finite difference evaluation, assuming that gas0 is valid state already.
   if not tol then tol = 0.0001 end
   local gas1 = GasState:new{gmodel}
   copyValues(gas0, gas1)
   local p0 = gas0.p
   local rho0 = gas0.rho
   local e0 = gas0.e[1]
   --
   local drho = rho0 * tol
   gas1.rho = rho0 + drho
   gmodel:updateThermoFromRHOE(gas1)
   local dpdrho = (gas1.p - p0)/drho
   --
   gas1.rho = rho0
   local de = e0 * tol
   gas1.e[1] = e0 + de
   gmodel:updateThermoFromRHOE(gas1)
   local dpde = (gas1.p - p0)/de
   --
   return dpdrho, dpde
end

------------------------------------------------------------------
-- Start the main script...
debug = false
do_gas_dynamic_accommodation = true
print("# Reacting pipe flow -- Bittker-Scullin test case 3.")
print(sample_header)

-- Gas model setup
-- nsp, nmodes, gmodel = setGasModel("h2-o2-n2-9sp.lua")
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
gas0.T[1] = 1559.0 -- degree K
molefInit = {O2=0.1480, N2=0.5562, H2=0.2958}
massfInit = gmodel:molef2massf(molefInit)
for i=1, nsp do
   local name = gmodel:speciesName(i)
   gas0.massf[name] = massfInit[name]
end
gmodel:updateThermoFromPT(gas0)
local dt_suggest = 1.0e-8  -- suggested starting time-step for chemistry updater
print(sample_data(x, v, gas0, dt_suggest))

print("# Start reactions...")
chemUpdate = ChemistryUpdate:new{filename="h2-o2-n2-9sp-18r.lua", gasmodel=gmodel}
local t = 0 -- time is in seconds
local t_final = 22.0e-6
local t_inc = 0.05e-6
local nsteps = math.floor(t_final / t_inc)
for j=1, nsteps do
   -- At the start of the step...
   local rho = gas0.rho
   local T = gas0.T[1]
   local p = gas0.p
   local e = gas0.e[1]
   --
   -- Do the chemical increment.
   local gas1 = GasState:new{gmodel} -- make the new one as a clone
   copyValues(gas0, gas1)
   dt_suggest = chemUpdate(gas1, t_inc, dt_suggest, gmodel)
   gmodel:updateThermoFromRHOE(gas1)
   --
   local de_chem = gas1.e[1] - e
   local dp_chem = gas1.p - p
   if debug then print("# de_chem=", de_chem, "dp_chem=", dp_chem) end
   --
   if do_gas_dynamic_accommodation then
      -- Do the gas-dynamic accommodation after the chemical change.
      local Etot = e + 0.5*v*v
      local dfdr, dfde = eos_derivatives(gas1, gmodel)
      if debug then print("# dfdr=", dfdr, "dfde=", dfde) end
      --[[ The Python code used a linear solve to get the accommodation increments.
	 A = np.array([ [v,      rho,        0.0, 0.0  ],
                      [0.0,    rho*v,      1.0, 0.0  ],
                      [v*Etot, rho*Etot+p, 0.0, rho*v],
                      [dfdr,   0.0,       -1.0, dfde ] ]);
        b = np.array([0.0, -dp_chem, -rho*v*de_chem, 0.0])
        dq = np.linalg.solve(A,b)
        drho, dv, dp_gda, de_gda = dq
      --]]
      -- The Lua code computes the accommodation increments explicitly,
      -- using expressions from Maxima.
      local denom = rho*rho*v*v - dfdr*rho*rho - dfde*p
      local drho = (dp_chem - de_chem*dfde)*rho*rho / denom
      local dv = -(dp_chem - de_chem*dfde)*rho*v / denom
      local dp_gda = -(de_chem*dfde*rho*rho*v*v - dfdr*dp_chem*rho*rho 
			  - dfde*dp_chem*p) / denom
      local de_gda = -(de_chem*rho*rho*v*v - de_chem*dfdr*rho*rho - dp_chem*p) / denom
      if debug then 
	 print("# drho=", drho, "dv=", dv, "dp_gda=", dp_gda, "de_gda=", de_gda)
	 print("# residuals=", v*drho + rho*dv, rho*v*dv + dp_gda + dp_chem,
	       v*Etot*drho + (rho*Etot+p)*dv + rho*v*de_gda + rho*v*de_chem,
	       dfdr*drho - dp_gda + dfde*de_gda)
      end
      -- Add the accommodation increments.
      gas1.rho = gas0.rho + drho
      v1 = v + dv
      p1_check = gas1.p + dp_gda
      gas1.e[1] = gas1.e[1] + de_gda
      gmodel:updateThermoFromRHOE(gas1)
      if debug then
	 print("# At new point for step ", j, ": gas1.p=", gas1.p, "p1_check=", p1_check, 
	       "rel_error=", math.abs(gas1.p-p1_check)/p1_check)
      end
   else
      -- Don't do the gas-dynamic accommodation.
      v1 = v
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
