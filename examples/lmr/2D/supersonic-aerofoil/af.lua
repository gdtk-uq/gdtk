-- "Kindwind", diamond aerofoil with trailing flap, for lmr.
-- @author: NNG and BD, June 2024 

config.solver_mode = "steady"
config.dimensions = 2
config.flux_calculator = "ausmdv"
config.viscous = true
config.turbulence_model = "spalart_allmaras_edwards"
config.with_local_time_stepping=true

--config.unstructured_limiter = "svan_albada"
config.unstructured_limiter = "hvenkat"
config.extrema_clipping=false                     -- On for safety, off for deep convergence
config.smooth_limiter_coeff = 1e-2
config.suppress_reconstruction_at_boundaries = true

-- Gas model and flow conditions Mach 7 Enthalpy
nsp, nmodes, gm = setGasModel('lmrsim/gm-air.lua')
p_inf = 2624.0 -- Pa
u_inf = 2221.0 -- m/s
T_inf = 250.0  -- K

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf
gm:updateThermoFromPT(gas_inf); gm:updateSoundSpeed(gas_inf); gm:updateTransCoeffs(gas_inf)

-- Turbulence quantities estimate
turb_lam_viscosity_ratio = 5.0 -- Transitional starting ratio from LARC website
nu_inf = gas_inf.mu/gas_inf.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf

dofile('angle_of_attack.lua')
theta = math.rad(angle_in_degrees)

inflow = FlowState:new{p=p_inf, velx=u_inf*math.cos(theta), vely=u_inf*math.sin(theta), T=T_inf, nuhat=nuhat_inf}
flowDict = {}
flowDict["inflow"] = inflow

bcDict = {inflow=InFlowBC_Supersonic:new{flowCondition=inflow},
             upper=InFlowBC_Supersonic:new{flowCondition=inflow},
             lower=InFlowBC_Supersonic:new{flowCondition=inflow},
             wall=WallBC_NoSlip_FixedT:new{Twall=300.0, group="wall"},
             outflow=OutFlowBC_Simple:new{},
             METIS_INTERIOR=ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}

makeFluidBlocks(bcDict, flowDict)

-- loads settings
config.boundary_groups_for_loads = "wall"

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 1200,
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
   write_loads = true,
   total_snapshots = 2,
   steps_between_snapshots = 100,
   steps_between_diagnostics = 1,
   steps_between_loads_update = 100,
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   frozen_preconditioner = true,
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 25,
   linear_solve_tolerance = 0.01,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.1,
   start_cfl = 2.0,
   max_cfl = 1.0e4,
   auto_cfl_exponent = 0.8
}
