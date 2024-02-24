-- transient.lua
print("Set up transient solve of Mach 1.5 flow over a 20 degree cone.")
--
-- 0. Assume that a previous processing has step set up the grid.
--
-- 1. Domain type, gas model and flow states
config.solver_mode = "transient"
config.axisymmetric = true
setGasModel('ideal-air.gas')
initial = FlowState:new{p=5955.0, T=304.0} -- Pa, degrees K
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
flowDict = {initial=initial, inflow=inflow}
--
-- 2. Fluid blocks, with initial flow states and boundary conditions.
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{},
   wall=WallBC_WithSlip:new{group="wall"}
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- 4. Simulation parameters.
config.flux_calculator= "ausmdv"
config.interpolation_order = 2
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
-- config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5
config.extrema_clipping = false

config.boundary_groups_for_loads = "wall"
