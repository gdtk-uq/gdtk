--
-- Supersonic Flow Through a Duct with Triangular Bumps at Mach 2.0 taken from ref. [1].
--
-- references:
-- [1] New Unstructured-Grid Limiter Functions,
--     H. Nishikawa,
--     AIAA SciTech Forum, 2022
--
-- @author: Kyle A. Damm (2024-08-06)
--

-- ==========================================================
-- General settings
-- ==========================================================
config.solver_mode          = "steady"
config.field_format         = "rawbinary"
config.grid_format          = "rawbinary"

-- ==========================================================
-- Freestream conditions
-- ==========================================================
nsp, nmodes, gmodel = setGasModel('ideal-air.gas')
gs = GasState:new{gmodel}

-- define free stream condition from ref. [1].
M_inf  = 2.0
gs.rho = 1.0   -- kg/m^3
gs.T   = 300.0 -- K
gmodel:updateThermoFromRHOT(gs)
gmodel:updateSoundSpeed(gs)
V_inf = M_inf*gs.a

-- set inflow and initial condition
flowDict            = {}
inflow              = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}
flowDict["inflow"]  = inflow
initial             = inflow
flowDict["initial"] = initial

-- ==========================================================
-- Define the flow domain using an imported grid
-- ==========================================================
grid_location       = "./"
config.dimensions   = 2
nblocks             = 1

-- import grids
grids = {}
for i=0,nblocks-1 do
   fileName = string.format(grid_location.."block_%d_duct.su2", i)
   grids[i] = registerFluidGrid{
      grid=UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1},
      fsTag="initial"
   }
end

-- define boundary conditions
bcDict = {
   inflow         = InFlowBC_Supersonic:new{flowCondition=inflow},
   outflow        = OutFlowBC_SimpleExtrapolate:new{},
   wall           = WallBC_WithSlip:new{},
   METIS_INTERIOR = ExchangeBC_MappedCell:new{cell_mapping_from_file=true,
                                              filename=grid_location.."mapped_cells",
                                              list_mapped_cells=false}
}

-- make fluidblocks
makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver settings
-- ==========================================================

-- invsicid flux settings
config.flux_calculator       = "roe"
config.interpolation_order   = 2
config.extrema_clipping      = false
config.thermo_interpolator   = "rhop"
config.unstructured_limiter  = "venkat_mlp"
config.smooth_limiter_coeff  = 1.0
config.use_extended_stencil  = true

-- Set temporal integration settings
NewtonKrylovGlobalConfig{
   -- phases
   number_of_phases = 3,
   max_steps_in_initial_phases = { 500, 500 },
   phase_changes_at_relative_residual = { 5.0e-02, 1.0e-06 },

   -- Preconditioner
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 1,

   -- nonlinear system solver settings
   max_newton_steps = 100000,
   max_consecutive_bad_steps = 10,
   stop_on_relative_residual = 1.0e-10,

   -- linear system solver settings
   frechet_derivative_perturbation = 1.0e-50,
   use_scaling = true,
   max_linear_solver_iterations = 100,
   max_linear_solver_restarts = 0,

   -- continuation settings
   inviscid_cfl_only = true,
   use_line_search = true,
   line_search_order = 3,
   use_physicality_check = true,
   allowable_relative_mass_change = 0.5,
   min_relaxation_factor_for_update = 0.1,
   min_relaxation_factor_for_cfl_growth = 0.5,

   -- output settings
   number_of_steps_for_setting_reference_residuals = 0,
   steps_between_status = 1,
   total_snapshots = 5,
   steps_between_snapshots = 10,
   steps_between_diagnostics = 1,
   write_limiter_values = true,
   write_gradient_values = true,
   write_residual_values = true
}

NewtonKrylovPhase:new{
   -- preconditioner settings
   frozen_preconditioner = true,
   use_adaptive_preconditioner = true,
   steps_between_preconditioner_update = 5,

   -- linear system solver settings
   linear_solve_tolerance = 1.0e-01,
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 1,

   -- cfl settings
   use_auto_cfl = true,
   use_local_timestep = true,
   threshold_relative_residual_for_cfl_growth = 0.99,
   start_cfl = 5.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}

NewtonKrylovPhase:new{
   -- linear system solver settings
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = -1
}

NewtonKrylovPhase:new{
   frozen_limiter_for_residual = true,
   frozen_limiter_for_jacobian = true,
   start_cfl = -1
}
