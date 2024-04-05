print("EAST facility with nozzle and 2T nitrogen")
-- Peter J. 2024-04-05
--
nsp, nmodes, gm = setGasModel("two-temp-n2.gas")
config.reacting = true
config.reactions_file = "chem.chem"
config.energy_exchange_file = "VT-relaxation-time-selection.kin"
--
-- Conditions before shock and at throat.
initial = FlowState:new{p=20.0e3, T=300.0, T_modes=300.0}
inflow = FlowState:new{p=5.898e6, T=5135.2, T_modes=5135.2, velx=1382.8}
--
flowDict = {
   initial=initial,
   inflow=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
makeFluidBlocks(bcDict, flowDict)
--
setHistoryPoint{x=0.095, y=0.0}  -- corresponds to c0, nozzle exit
--
config.solver_mode = "transient"
config.max_time = 225.0e-6 -- s
config.max_step = 30000
config.cfl_value = 0.4
config.dt_init = 1.0e-9
config.dt_plot = 10.0e-6
config.dt_history = 0.1e-6
