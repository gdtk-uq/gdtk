print("Under-expanded jet of nitrogen with pressure ratio 40.")
-- PJ 2021-03-27 Eilmer4 example
--    2024-03-27 Eilmer5 port
--
config.dimensions = 2
config.axisymmetric = true
config.solver_mode = "transient"

-- Gas model and flow conditions.
setGasModel('ideal-n2.gas')
tank = FlowState:new{p=101.325e3, T=300.0}
sspeed = tank.a
print("sound speed=", sspeed)
reservoir = FlowState:new{p=40.0*tank.p, T=300.0}
--
flowDict = {
   reservoir = reservoir,
   tank = tank
}
bcDict = {
   ambient = InOutFlowBC_Ambient:new{flowState=tank},
   outflow = OutFlowBC_Simple:new{}
}
makeFluidBlocks(bcDict, flowDict)
--
-- History points for reservoir pressure and pitot probe.
setHistoryPoint{x=-0.020, y=0.0275}
setHistoryPoint{x=0.060, y=0.010}
config.dt_history = 10.0e-6
--
mpiDistributeBlocks{ntasks=7} -- ntasks needs to match the mpirun command.
config.max_time = 5.0e-3  -- seconds
config.max_step = 100000
config.dt_plot = 0.050e-3
