--
-- Conjuagte heat transfer simulation of hypersonic flow over a hollow cylinder with finite thickness,
-- freestream condition and experimental data digitized from ref. [1].
--
-- references:
-- [1] Experimental Study of Shock Wave Interference Heating on a Cylindrical Leading Edge
--     A. R. Wieting
--     NASA Technical Memorandum, TM-100484
--
-- author: Kyle A. Damm
-- date:   2023-11-09
--

-- ==========================================================
-- General settings
-- ==========================================================
job_title = "Mach 6.47 air flow over a hollow cylinder."
print(job_title)
config.title       = job_title
config.print_count = 1
config.flow_format = "eilmer4binary"
config.new_flow_format      = true
config.save_residual_values = true

-- ==========================================================
-- Freestream conditions
-- ==========================================================
-- select gas model
nsp, nmodes, gmodel = setGasModel('ideal-air.gas')
gs = GasState:new{gmodel}

-- run 37 conditions taken taken from Appendix D in ref. [1]
Mach       = 6.47
T_wall     = 294.444   -- K
gs.T       = 241.5     -- K
gs.p       = 648.11    -- Pa
V_inf      = 2034.235  -- m/s

-- update some gas properties
gmodel:updateThermoFromPT(gs)
gmodel:updateTransCoeffs(gs)
gmodel:updateSoundSpeed(gs)

-- set inflow condition
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}

-- turbulence model settings
config.turbulence_model  = "none"
if config.turbulence_model == "none" then
   -- do nothing (laminar flow doesn't need any turbulence variables set).
elseif config.turbulence_model == "spalart_allmaras" then
   turb_lam_viscosity_ratio = 5.0
   nu_inf                   = gs.mu/gs.rho
   nuhat_inf                = turb_lam_viscosity_ratio*nu_inf
   inflow.turb[1]           = nuhat_inf
elseif config.turbulence_model == "k_omega" or config.turbulence_model == "k_omega_vorticity" then
   turb_intensity             = 0.1
   turb_to_laminar_visc_ratio = 1.0
   tke_inf                    = (1.2/2.0) * (turb_intensity/100 * V_inf)^2
   nu_inf                     = gs.mu/gs.rho
   omega_inf                  = (tke_inf/nu_inf)/turb_to_laminar_visc_ratio
   inflow.turb[1]             = tke_inf
   inflow.turb[2]             = omega_inf
else
   print("WARNING: undefined behaviour for turbulence model.")
end

-- set initial condition
initial = inflow

-- ==========================================================
-- Define the fluid and solid domain using an imported grid
-- ==========================================================
config.coupling_with_solid_domains = "steady_fluid_transient_solid"
config.dimensions   = 2
config.axisymmetric = false
-- Set up the geometry for defining the grid
Ro = 3.81e-02 -- outer nose radis, metres
Ri = 2.52e-02 -- inner nose radius, metres

-- define FluidBlocks
n0i = 100
n0j = 50

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-Ro, y=0.0}
c = Vector3:new{x=0.0, y=Ro}
d = { Vector3:new{x=-1.5*Ro, y=0.0}, Vector3:new{x=-1.5*Ro, y=Ro},
      Vector3:new{x=-Ro, y=2*Ro}, Vector3:new{x=0.0, y=3*Ro} }

sphere_edge = Arc:new{p0=b, p1=c, centre=a}
psurf0 = makePatch{north=Line:new{p0=d[#d], p1=c}, south=Line:new{p0=d[1], p1=b},
		  east=sphere_edge, west=Bezier:new{points=d}}
cf_radial0 = RobertsFunction:new{end0=false, end1=true, beta=1.2}

grid0 = StructuredGrid:new{psurface=psurf0, niv=n0i+1, njv=n0j+1,
			  cfList={north=cf_radial0, south=cf_radial0}}

blk0 = FBArray:new{grid=grid0, initialState=initial, label='fluidblk',
                       bcList={north=OutFlowBC_SimpleExtrapolate:new{},
                               east=WallBC_AdjacentToSolid2:new{},
                               south=WallBC_WithSlip:new{},
                               west=InFlowBC_Supersonic:new{flowState=inflow}},
                       nib=2, njb=2}

-- define SolidBlocks
n1i = 25

e = Vector3:new{x=-Ri, y=0.0}
f = Vector3:new{x=0.0, y=Ri}

sphere_outer_edge = Arc:new{p0=b, p1=c, centre=a}
sphere_inner_edge = Arc:new{p0=e, p1=f, centre=a}
psurf1 = makePatch{north=Line:new{p0=c, p1=f}, south=Line:new{p0=b, p1=e},
		  east=sphere_inner_edge, west=sphere_outer_edge}
cf_radial1 = RobertsFunction:new{end0=true, end1=false, beta=1.1}

grid1 = StructuredGrid:new{psurface=psurf1, niv=n1i+1, njv=n0j+1,
                           cfList={north=cf_radial1, south=cf_radial1}}

blk1 = SolidBlockArray{grid=grid1, initTemperature=T_wall, label='solidblk',
                       properties={rho=8030.0, k=16.24, Cp=502.48},
                       bcList={north=SolidAdiabaticBC:new{},
                               east=SolidFixedTBC:new{T=Twall},
                               south=SolidConstantFluxBC:new{},
                               west=SolidAdjacentToGasBC2:new{}},
                       nib=2, njb=2}

identifyBlockConnections()

-- every MPI task needs a FluidBlock, so make sure you distribute the blocks accordingly
-- e.g. here we have 4 FluidBlocks + 4 SolidBlocks divided amongst 4 MPI tasks
mpiTasks = mpiDistributeBlocks{ntasks=4, dist="load-balance", preassign={[0]=1}}

-- ==========================================================
-- Solver settings
-- ==========================================================

-- loads settings (adjacent_to_solid is a special internal boundary name for coupled fluid/solid boundaries)
config.boundary_groups_for_loads = "adjacent_to_solid"
config.write_loads = true

-- solid domain settings
config.solid_domain_augmented_deriv_avg = false
config.gasdynamic_update_scheme         = 'rkl1'
config.cfl_value = 0.5
config.dt_init   = 1.0e-12
config.max_time  = 1.0
config.max_step  = 1e8
config.dt_plot   = config.max_time

-- invsicid flux settings
config.flux_calculator       = "adaptive_hanel_ausmdv"
config.apply_entropy_fix     = false
config.strict_shock_detector = true
config.compression_tolerance = -0.05
config.shear_tolerance       = 0.01
config.interpolation_order   = 2
config.thermo_interpolator   = "rhop"
config.extrema_clipping      = false
config.apply_limiter         = true
config.apply_heuristic_pressure_based_limiting = true

-- viscous flux settings
config.viscous            = true
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
-- special blending to improve start-up stability for k-omega model
if config.turbulence_model == "k_omega" or config.turbulence_model == "k_omega_vorticity" then
   config.diffuse_wall_bcs_on_init = true
   config.number_init_passes = 10
else
   config.diffuse_wall_bcs_on_init = false
end

config.with_local_time_stepping = true
SteadyStateSolver{
   -- preconditioner settings
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 1,
   start_preconditioning = 0,
   preconditioner_sigma = 1.0e-50,
   use_adaptive_preconditioner = true,

   -- linear solve general settings
   use_scaling = true,
   use_physicality_check = true,
   physicality_check_theta = 0.25,
   use_complex_matvec_eval = true,

   -- freeze limiter for unstructured grids (and shock detector for all grids)
   limiter_freezing_residual_reduction = 1.0e-05,
   limiter_freezing_count = 1,

   -- stopping criteria
   cfl_min = 1.0e-04,
   number_pre_steps = 0,
   number_total_steps = 1000,
   stop_on_relative_global_residual = 1.0e-08,

   -- settings for GMRES iterative solver
   max_outer_iterations = 50,
   max_restarts = 0,

   -- CFL settings
   include_turb_quantities_in_residual = false,
   residual_based_cfl_scheduling = true,
   inviscid_cfl = true,
   cfl_max = 1e6,

   -- start-up phase settings
   number_start_up_steps = 0,
   sigma0 = 1.0e-50,
   cfl0 = 1.0,

   -- final phase settings
   sigma1 = 1.0e-50,
   cfl1 = 1.0,
   p1 = 0.75,
   tau1 = 1.0,
   eta1 = 1.0e-02,
   eta_strategy = "constant",

   -- output settings
   snapshots_count = 1,
   number_total_snapshots = 10,
   write_diagnostics_count = 1,
   write_loads_count = 50,
}
