-- ramp15.lua
-- Reacting flow over a 15 degree ramp
-- Kyle Damm & Nick Gibbons
-- 2024-04-08

config.flux_calculator = "ausmdv"
config.interpolation_order = 1

nsp, nmodes, gmodel = setGasModel('lmrsim/species.lua')
config.reacting = true
config.reactions_file = 'lmrsim/reactions.lua'
config.solver_mode = 'steady'
--config.solver_mode = 'transient'

-- inflow gas state
Mach = 3.4
gs = GasState:new{gmodel}
gs.T = 700.0 -- Kelvin
gs.p = 101e3 -- Pa, 1 atm
gs.massf = {H2=0.005, O2=0.24, N2=0.755} -- general approximation of air
gamma = gmodel:gamma(gs)
R = gmodel:R(gs)
vel_inf = Mach*math.sqrt(gs.T*gamma*R)
inflowState = FlowState:new{p=gs.p, T=gs.T, velx=vel_inf, massf=gs.massf}

-- Define the flow domain using an imported grid.
--nblocks=8
--grids = {}
--for i=0,nblocks-1 do
--   fileName = string.format("su2-grid/block_%d_ramp15.su2", i)
--   registerFluidGrid{grid=UnstructuredGrid:new{filename=fileName, fmt="su2text"}, fsTag="inflow"}
--end

bcDict = {north = OutFlowBC_Simple:new{},
          east  = OutFlowBC_Simple:new{},
          south = WallBC_WithSlip:new{},
          west  = InFlowBC_Supersonic:new{flowCondition=inflowState},
          METIS_INTERIOR   = ExchangeBC_MappedCell:new{cell_mapping_from_file=true,
                                                       list_mapped_cells=false}}

flowDict = {inflow = inflowState}

makeFluidBlocks(bcDict, flowDict)

-- Transient timing settings
flow_length = 1.5
flow_time = flow_length/vel_inf
no_flow_times = 4.0

config.gasdynamic_update_scheme = "pc"
config.max_time = no_flow_times*flow_time
config.max_step = 80000
config.dt_init = 1.0e-10
config.cfl_value = 0.50
config.dt_plot = config.max_time/10.0

-- Steady timing settings
NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 2000,
   stop_on_relative_residual = 1.0e-12,
   number_of_phases = 1,
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = true,
   max_linear_solver_iterations = 80,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-30,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-30,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 2,
   steps_between_snapshots = 100,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   frozen_preconditioner = true,
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 20,
   linear_solve_tolerance = 0.001,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.1,
   start_cfl = 1.0,
   max_cfl = 1.0e2,
   auto_cfl_exponent = 0.6
}

