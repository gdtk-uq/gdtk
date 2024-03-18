-- steady.lua
print("Set up steady-state solve of Mach 1.5 flow over a 20 degree cone.")
--
-- 0. Assume that a previous processing step has set up the grid.
--
-- 1. Domain type, gas model and flow states
config.solver_mode = "steady"
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
   wall=WallBC_WithSlip:new{group="wall"},
   mapped_cells = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, symmetric_mapping=true}
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- 3. Simulation parameters.
config.flux_calculator= "ausmdv"
config.interpolation_order = 2

config.boundary_groups_for_loads = "wall"

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-6,
   number_of_phases = 2,
   max_steps_in_initial_phases = { 15 },
   use_physicality_check = true,
   max_linear_solver_iterations = 10,
   total_snapshots = 3,
   steps_between_status = 1,
   steps_between_snapshots = 5,
   steps_between_diagnostics = 1,
   write_loads = true
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.9,
   start_cfl = 2.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 0.9
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 10.0
}
