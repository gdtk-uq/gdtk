-- poshax.lua
-- Post-shock-relaxation calculation using the general gas model.
--
-- PJ, 2018-04-21
-- Started with Rowan's description in Section 4.4 of his PhD thesis
-- and then built upon the gas-model API, the kinetics API,
-- the quasi-1D gas flow functions, and linearized 1D flow constraints
-- as used in the duct-with-hydrogen-combustion example.
-- This is a proof of concept code.
--
-- Typical use:
-- $ e4shared --custom-post --script-file=poshax.lua
--
-------------------------------------------------------------------------
print("Initialise a gas model.")
print("air 7 species Gupta-et-al reactions")
gmodel = GasModel:new{"air-7sp-gas-model.lua"}
chemUpdate = ChemistryUpdate:new{filename="air-7sp-chemistry.lua",
                                 gasmodel=gmodel}
-- The example here matches the case discussed on page 63 of the thesis.
state1 = GasState:new{gmodel}
state1.p = 133.3 -- Pa
state1.T = 300.0 -- degree K
state1.massf = {N2=0.78, O2=0.22}
print("Free stream conditions, before the shock.")
gmodel:updateThermoFromPT(state1)
gmodel:updateSoundSpeed(state1)
print("state1:"); printValues(state1)
mach1 = 12.28
V1 = mach1 * state1.a
print("mach1:", mach1, "V1:", V1)

print("Stationary normal shock with chemically-frozen gas.")
state2, V2, Vg = gasflow.normal_shock(state1, V1)
print("    V2=", V2, "Vg=", Vg)
print("    state2:"); printValues(state2)

-----------------------------------------------------------------------
print("Decide what to write to the output file.")
output_filename = "air-Mach_12.3.data"
write_massfractions = true
write_molefractions = true

speciesNames = {}
for i=1,gmodel:nSpecies() do
   -- remember that Dlang starts counting from zero
   speciesNames[i] = gmodel:speciesName(i-1)
end

-- Output data format is a little like the old poshax code.
header = "# 1:x(m) 2:T(degK) 3:p(Pa) 4:rho(kg/m**3) 5:e(J/kg) 6:v(m/s)"
columnNum = 7
for i,name in ipairs(speciesNames) do
   if write_massfractions then
      header = header .. string.format(" %d:massf_%s", columnNum, name)
      columnNum = columnNum + 1
   end
   if write_molefractions then
      header = header .. string.format(" %d:molef_%s", columnNum, name)
      columnNum = columnNum + 1
   end
end
header = header .. string.format(" %d:dt_suggest(s)\n", columnNum)

function string_data(x, v, gas, dt_suggest)
   local str = string.format("%g %g %g %g %g %g",
                             x, gas.T, gas.p, gas.rho, gas.u, v)
   local conc = gmodel:massf2conc(gas)
   for i,name in ipairs(speciesNames) do
      if write_massfractions then
         str = str .. string.format(" %g", gas.massf[name])
      end
      if write_molefractions then
         str = str .. string.format(" %g", conc[name])
      end
   end
   str = str .. string.format(" %g\n", dt_suggest)
   return str
end

---------------------------------------------------------------------------
print("Relaxing flow starts here.")
debug = false

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
--
print("Step in time, allowing the gas to react and drift along in x.")
gas0 = GasState:new{gmodel}
copyValues(state2, gas0) -- start with post-shock conditions
x = 0 -- position of shock, m
x_final = 1.0 -- the marching will stop at this point, metres
v = V2 -- velocity of gas post-shock, m/s
dt_suggest = 1.0e-8  -- suggested starting time-step for chemistry updater
--
f = io.open(output_filename, 'w')
f:write(header)
f:write(string_data(x, v, gas0, dt_suggest))
--
t = 0 -- time of drift is in seconds
t_final = 3.0e-3 -- User-selected drift time, s
t_inc = 0.1e-6 -- User-selected time steps, s
nsteps = math.floor(t_final / t_inc)
-- Notes:
--  Select t_final to be long enough to reach x_final.
--  Select t_inc small enough to capture fast changes near x=0.
--
for j=1, nsteps do
   -- At the start of the step, we have GasState 0.
   -- Make some shorter names.
   local rho, T, p, u = gas0.rho, gas0.T, gas0.p, gas0.u
   --
   -- Do the chemical increment.
   -- Make the new GasState as a clone and then update it.
   local gas1 = GasState:new{gmodel}
   copyValues(gas0, gas1)
   dt_suggest = chemUpdate:updateState(gas1, t_inc, dt_suggest, gmodel)
   gmodel:updateThermoFromRHOU(gas1)
   --
   local du_chem = gas1.u - u
   local dp_chem = gas1.p - p
   if debug then print("du_chem=", du_chem, "dp_chem=", dp_chem) end
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
   --]=]
   -- Compute the accommodation increments using expressions from Maxima.
   local denom = rho*rho*v*v - dfdr*rho*rho - dfdu*p
   local drho = (dp_chem - du_chem*dfdu)*rho*rho / denom
   local dv = -(dp_chem - du_chem*dfdu)*rho*v / denom
   local dp_gda = -(du_chem*dfdu*rho*rho*v*v - dfdr*dp_chem*rho*rho 
                       - dfdu*dp_chem*p) / denom
   local du_gda = -(du_chem*rho*rho*v*v - du_chem*dfdr*rho*rho - dp_chem*p) / denom
   if debug then 
      print("drho=", drho, "dv=", dv, "dp_gda=", dp_gda, "du_gda=", de_gda)
      print("residuals=", v*drho + rho*dv, rho*v*dv + dp_gda + dp_chem,
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
      print("At new point for step ", j, ": gas1.p=", gas1.p, "p1_check=", p1_check, 
            "rel_error=", math.abs(gas1.p-p1_check)/p1_check)
   end
   -- Have now finished the chemical and gas-dynamic update.
   t = t + t_inc
   x = x + 0.5*(v + v1) * t_inc
   f:write(string_data(x, v1, gas1, dt_suggest))
   if x > x_final then break end
   -- House-keeping for the next step.
   v = v1
   copyValues(gas1, gas0) -- gas0 will be used in the next iteration
end
f:close()
print("Done stepping.")

