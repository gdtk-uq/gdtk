print("Prepare FluidBlocks and config for forward-facing-step.")
-- PJ 2021-04-21 Eilmer4 staged prep
--    2024-03-26 port to Eilmer 5
--
config.dimensions = 2
config.solver_mode = "transient"
config.max_time = 5.0e-3  -- seconds
config.max_step = 6000
config.dt_plot = 1.0e-3
--
setGasModel('ideal-air.gas')
initial = FlowState:new{p=101.325e3, T=300.0}
print("sound speed=", initial.a)
inflow = FlowState:new{p=101.325e3, T=300.0, velx=3.0*initial.a, vely=0.0}
--
flowDict = {
   INFLOW = inflow
}
bcDict = {
   INFLOW = InFlowBC_Supersonic:new{flowState=inflow},
   OUTFLOW = OutFlowBC_Simple:new{}
}
makeFluidBlocks(bcDict, flowDict)
--
mpiDistributeBlocks{ntasks=3} -- ntasks needs to match the mpirun command.
