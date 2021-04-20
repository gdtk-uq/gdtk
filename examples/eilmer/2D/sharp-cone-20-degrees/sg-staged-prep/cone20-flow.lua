-- cone20-flow.lua
-- Simulation specification, following a grid specification.
-- 2021-04-19  PeterJ
--
config.title = "Mach 1.5 flow over a 20 degree cone."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_plot = 1.5e-3
--
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
flowDict = {}
flowDict["initial-gas"] = FlowState:new{p=5955.0, T=304.0}
flowDict["inflow-gas"] = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
--
bcDict = {
   INFLOW=InFlowBC_Supersonic:new{flowState=inflow},
   OUTFLOW=OutFlowBC_Simple:new{},
}
--
makeFluidBlocks(bcDict, flowDict)
