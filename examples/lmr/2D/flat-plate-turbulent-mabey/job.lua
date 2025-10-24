--
-- Turbulent flow over a flat plate, details taken from ref. [1],
-- adapted from the Eilmer4 script.
--
-- references:
--
-- [1] Mabey test case (AGARDograph 223 - Test series 7402),
--     (Referenced from Fernholz & Finley (1977),
--     AGARDograph No. 223, "A critical compilation of
--     compressible turbulent boundary layer data.")
--
-- author: Kyle A. Damm
-- date:   2025-10-01
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
nsp, nmodes, gm = setGasModel('ideal-air.gas')
gs = GasState:new{gm}

-- define flow state from ref. [1] to match Mabey's data set 74021802
gs.T       = 61.16     -- degrees K
gs.p       = 3.16e3    -- Pa
V_inf      = 712.9     -- m/s

-- iupdate thermodynamic transfer coefficients
gm:updateThermoFromPT(gs)
gm:updateSoundSpeed(gs)
gm:updateTransCoeffs(gs)

-- Estimate turbulence quantities for free stream by specifying
-- turbulence intensity and turbulent-to-laminar viscosity ratio
config.turbulence_model = "k_log_omega"
turb_intensity = 0.01
turb_lam_viscosity_ratio = 1.0
tke_inf = 1.5 * (turb_intensity * V_inf)^2
mu_t_inf = turb_lam_viscosity_ratio * gs.mu
omega_inf = gs.rho * tke_inf / mu_t_inf
if string.find(config.turbulence_model, "log") then
   omega_inf = math.log(omega_inf)
end

-- set inflow and initial conditions
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V_inf,  tke=tke_inf, omega=omega_inf}
initial = inflow
print("Inflow Check: \n", inflow)

-- ==========================================================
-- Define the flow domain using an imported grid
-- ==========================================================
grid_location       = "./su2grid/"
config.dimensions   = 2
config.axisymmetric = false
nblocks             = 4

-- import grids
grids = {}
for i=0,nblocks-1 do
   fileName = string.format(grid_location.."block_%d_grid.su2", i)
   grids[i] = registerFluidGrid{
      grid=UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1},
      fsTag="initial"
   }
end

-- define boundary conditions
bcDict = {
   inflow         = InFlowBC_Supersonic:new{flowCondition=inflow},
   outflow        = OutFlowBC_SimpleExtrapolate:new{},
   wall           = WallBC_NoSlip_Adiabatic0:new{group="wall"},
   METIS_INTERIOR = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, filename=grid_location.."mapped_cells"}
}

-- make fluidblocks
flowDict            = {}
flowDict["initial"] = initial
flowDict["inflow"]  = inflow
makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver settings
-- ==========================================================

-- loads settings
config.write_loads               = true
config.boundary_groups_for_loads = "wall"

-- invsicid flux settings
config.flux_calculator             = "ausmdv"
config.apply_entropy_fix           = false
config.interpolation_order         = 2
config.extrema_clipping            = false
config.thermo_interpolator         = "rhop"
config.unstructured_limiter        = "hvenkat"
config.smooth_limiter_coeff        = 1.0
config.inviscid_least_squares_type = "unweighted_qr"

-- viscous flux settings
config.viscous                    = true
config.spatial_deriv_calc         = "least_squares"
config.spatial_deriv_locn         = "cells"
config.viscous_least_squares_type = "weighted_qr"
config.include_boundary_faces_in_spatial_deriv_correction = true

-- Set temporal integration settings
NewtonKrylovGlobalConfig{
   -- phases
   number_of_phases = 3,
   max_steps_in_initial_phases = { 500, 500 },
   phase_changes_at_relative_residual = { 1.0e-03, 1.0e-6 },

   -- preconditioner settings
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,

   -- nonlinear system solver settings
   max_newton_steps = 2000,
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
   allowable_relative_mass_change = 0.9,
   min_relaxation_factor_for_update = 0.1,
   min_relaxation_factor_for_cfl_growth = 0.5,

   -- output settings
   number_of_steps_for_setting_reference_residuals = 5,
   steps_between_status = 1,
   --write_limiter_values = true,
   write_loads = true,
   total_snapshots = 5,
   steps_between_snapshots = 50,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   -- preconditioner settings
   frozen_preconditioner = true,
   use_adaptive_preconditioner = true,
   steps_between_preconditioner_update = 5,

   -- linear system solver settings
   linear_solve_tolerance = 1.0e-02,
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 1,
   use_residual_smoothing = true,

   -- cfl settings
   use_auto_cfl = true,
   use_local_timestep = true,
   threshold_relative_residual_for_cfl_growth = 0.9,
   start_cfl = 1.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 0.75
}

NewtonKrylovPhase:new{
   -- linear system solver settings
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   use_residual_smoothing = false,
   start_cfl = -1.0
}

NewtonKrylovPhase:new{
   frozen_shock_detector = true,
   frozen_limiter_for_residual = true,
   frozen_limiter_for_jacobian = true,
   start_cfl = -1
}
