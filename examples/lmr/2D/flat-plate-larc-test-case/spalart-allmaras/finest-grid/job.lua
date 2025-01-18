-- 2D Zero Pressure Gradient High Mach Number Flat Plate Validation Case
-- https://turbmodels.larc.nasa.gov/ZPGflatplateSS_val_sa.html
-- @author: Nick N. Gibbons (n.gibbons@uq.edu.au)
--
-- 2025-01-18, RJG
--    Ported to Eilmer 5
--

title = "LARC Mach 5.0 flow over a flat plate (spalart-allmaras): SIM PREPARATION"
print(title)

config.dimensions = 2
config.solver_mode = 'steady'
config.turbulence_model = "spalart_allmaras"
config.viscous = true
config.flux_calculator = "ausmdv"
config.interpolation_order = 2

config.unstructured_limiter = "svan_albada"
config.smooth_limiter_coeff = 1e-6
config.extrema_clipping = false
config.suppress_reconstruction_at_boundaries = true

config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.diffuse_wall_bcs_on_init = true
config.number_init_passes = 25

-- Gas model and flow conditions to match ZPG case M=5, Tw/Tinf=5.450
-- Shout out to this incredibly inconvenient problem description...
nsp, nmodes, gm = setGasModel('ideal-air.gas')
Re_unit = 15e6
T_inf = 540.6*5/9 -- (The 540.6 is degrees Rankine, convert to Kelvin)
M_inf = 5.0
Tw_on_Tinf = 5.450
Tw = Tw_on_Tinf*T_inf

-- Compute the gas state in physical units from its nondimensional description
gas_inf = GasState:new{gm}
gas_inf.T = T_inf
gm:updateSoundSpeed(gas_inf)
u_inf = M_inf*gas_inf.a
gm:updateTransCoeffs(gas_inf)
gas_inf.rho = Re_unit*gas_inf.mu/u_inf
gm:updateThermoFromRHOT(gas_inf)
p_inf = gas_inf.p

-- Use updated gas properties to estimate turbulence quantities
turb_lam_viscosity_ratio = 5.0 -- Fully turbulent, equation (8) from Allmaras (2012)
nu_inf = gas_inf.mu/gas_inf.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, nuhat=nuhat_inf}
--inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
print("Inflow Check\n", inflow)
print("gas_inf.rho: ", gas_inf.rho, "gas_inf.mu: ", gas_inf.mu)

bcDict = {
   supersonic=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{},
   fixedT=WallBC_NoSlip_FixedT0:new{Twall=Tw, wall_function=false, group="wall"},
   symm=WallBC_WithSlip:new{},
   mapped_cells = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, symmetric_mapping=true}
}
flowDict = {
   initial=inflow
}
makeFluidBlocks(bcDict, flowDict)

-- loads settings
config.boundary_groups_for_loads = "wall"
config.write_loads = true

config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.diffuse_wall_bcs_on_init = true
config.number_init_passes = 5
config.extrema_clipping = false


NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 500,
   stop_on_relative_residual = 1.0e-8,
   number_of_phases = 1,
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = true,
   max_linear_solver_iterations = 40,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 2,
   steps_between_snapshots = 100,
   steps_between_diagnostics = 1
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
   use_local_timestep = true,
   threshold_relative_residual_for_cfl_growth = 0.1,
   start_cfl = 2.0,
   max_cfl = 1.0e4,
   auto_cfl_exponent = 0.8
}

