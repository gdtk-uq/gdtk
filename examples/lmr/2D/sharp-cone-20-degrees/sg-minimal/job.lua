config.title = "Mach 1.5 flow over a 20 degree cone."
print(config.title)
--
-- 1. Domain type, gas model and flow states
config.axisymmetric = true
setGasModel('ideal-air.gas')
initial = FlowState:new{p=5955.0, T=304.0} -- Pa, degrees K
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
flowDict = {initial=initial, inflow=inflow}
--
-- 2. Geometry and grids
a0 = {x=0.0, y=0.0};     a1 = {x=0.0, y=1.0}
b0 = {x=0.2, y=0.0};     b1 = {x=0.2, y=1.0}
c0 = {x=1.0, y=0.29118}; c1 = {x=1.0, y=1.0}
--
quad0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
quad1 = AOPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
--
grid0 = registerGrid{
   grid=StructuredGrid:new{psurface=quad0, niv=11, njv=41},
   fsTag="inflow",
   bcTags={west="inflow"}
}
grid1 = registerGrid{
   grid=StructuredGrid:new{psurface=quad1, niv=31, njv=41},
   fsTag="initial",
   bcTags={east="outflow"}
}
identifyGridConnections()
--
-- 3. Fluid blocks, with initial flow states and boundary conditions.
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
--
makeFluidBlocks(bcDict, flowDict)
-- 4. Simulation parameters.
config.flux_calculator= "ausmdv"
config.interpolation_order = 2

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-6,
   number_of_phases = 2,
   phase_changes_at_steps = { 10 },
   use_physicality_check = true,
   max_linear_solver_iterations = 10,
   total_snapshots = 3,
   steps_between_status = 1,
   steps_between_snapshots = 5,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.9,
   start_cfl = 2.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 10.0
}
