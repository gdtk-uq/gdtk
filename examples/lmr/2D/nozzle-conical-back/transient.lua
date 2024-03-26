print("Flow through a conical nozzle.")
-- Conical nozzle from Back, Massier and Gier (1965)
-- Peter J. 2015-10-21 adpated from the Python version.
--          2024-03-26 ported to Eilmer 5.
--
setGasModel('ideal-air.gas')
fixed_pressure = false
if fixed_pressure then
   -- Directly specify the stagnation conditions for the subsonic inflow.
   stagnation_gas = FlowState:new{p=500.0e3, T=300.0}
   inBC = InFlowBC_FromStagnation:new{stagnationState=stagnation_gas,
                                      label="inflow-boundary"}
else
   -- We'll specify a inflow mass_flux (kg/s/m^^2) across the inlet
   -- and let the pressure be adjusted.
   -- Note that still need a stagnation condition and that the temperature is fixed.
   -- Note, also, that we start with a known wrong stagnation pressure.
   stagnation_gas = FlowState:new{p=400.0e3, T=300.0}
   inBC = InFlowBC_FromStagnation:new{stagnationState=stagnation_gas,
                                      mass_flux=275.16, relax_factor=0.2,
                                      label="inflow-boundary"}
end
low_pressure_gas = FlowState:new{p=30.0, T=300.0}

flowDict = {
   initial=low_pressure_gas,
   inflow=stagnation_gas
}
bcDict = {
   inflow=inBC,
   outflow=OutFlowBC_Simple:new{label="outflow-boundary"}
}
makeFluidBlocks(bcDict, flowDict)

config.solver_mode = "transient"
config.axisymmetric = true
config.flux_calculator = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"
config.max_time = 4.0e-3  -- seconds
config.max_step = 50000
config.dt_init = 1.0e-7
config.dt_plot = 0.2e-3

-- History locations near throat and exit in the second block
-- that models the supersonic part of the nozzle.
setHistoryPoint{ib=1, i=0, j=0}
setHistoryPoint{ib=1, i=9999, j=0} -- A large number gets the last cell in that direction.
config.dt_history = 10.0e-6
