--
-- This script is used to setup a simlulation of Rose and Stark's (ref. [1] & ref. [2])
-- hemispherical heat-transfer probe set in a shock tube.
--
-- references:
-- [1] Stagnation point heat-transfer measurements in dissociated air
--     P. H. Rose and W. I. Stark,
--     AIAA Journal (1957)
-- [1] Laminar heat transfer around blunt bodies in dissociated air
--     N. H. Kemp, P. H/ Rose, and R. W. Detra
--     AIAA Journal (1958)
--
-- @author: Kyle A. Damm (2024-08-23)
--          PAJ and RJG  (2016-08-28)
--

-- ==========================================================
-- General settings
-- ==========================================================
config.solver_mode          = "steady"

-- ==========================================================
-- Freestream conditions
-- ==========================================================
nsp, nmodes, gmodel = setGasModel('ideal-air.gas')
gs = GasState:new{gmodel}

-- define free stream condition for a M = 8 incident shock in air at 296 K (assuming fully-equilibrium chemistry)
gs.p  = 535.6  -- Pa
gs.T  = 2573.5 -- K
V_inf = 2436.5 -- m/s

-- set inflow and initial condition
flowDict            = {}
inflow              = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}
flowDict["inflow"]  = inflow
initial             = inflow
flowDict["initial"] = initial

-- ==========================================================
-- Define the flow domain using a native grid
-- ==========================================================
config.dimensions   = 2
config.axisymmetric = true
R = 6.6e-03 -- sphere radies in meters

-- define points
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-R, y=0.0}
c = Vector3:new{x=0.0, y=R}
d = { Vector3:new{x=-1.5*R, y=0.0}, Vector3:new{x=-1.5*R, y=R},
      Vector3:new{x=-R, y=2*R}, Vector3:new{x=0.0, y=3*R} }

-- define paths
p1 = Line:new{p0=d[#d], p1=c}
p2 = Arc:new{p0=b, p1=c, centre=a}
p3 = Line:new{p0=d[1], p1=b}
p4 = Bezier:new{points=d}

-- define surface
psurf = makePatch{north=p1, east=p2, south=p3, west=p4
}

-- define grid clustering
cluster = RobertsFunction:new{end0=false,end1=true,beta=1.05}

-- define grid
sgrid = StructuredGrid:new{psurface=psurf,
                           cfList = {north=cluster, east=None,
                                     south=cluster, west=None},
                           niv=200, njv=50}

grid = registerFluidGridArray{
   grid=sgrid,
   nib=1, njb=1,
   fsTag="initial",
   bcTags={north="outflow", east="wall", south="slip_wall", west="inflow"}
}

identifyGridConnections()

-- define boundary conditions
bcDict = {inflow         = InFlowBC_Supersonic:new{flowCondition=inflow},
          outflow        = OutFlowBC_SimpleExtrapolate:new{},
          wall           = WallBC_NoSlip_FixedT:new{Twall=296.0, group="wall"},
          mapped_cells=ExchangeBC_MappedCell:new{cell_mapping_from_file=true, symmetric_mapping=true}}

-- define fluidblocks
makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver settings
-- ==========================================================

-- loads settings
config.boundary_groups_for_loads = "wall"

-- invsicid flux settings
config.flux_calculator       = "ausmdv"
config.interpolation_order   = 2
config.extrema_clipping      = false
config.thermo_interpolator   = "rhop"
config.unstructured_limiter  = "hvenkat"
config.smooth_limiter_coeff  = 5.0
--config.apply_unstructured_limiter_stagnation_point_filter = true
--config.apply_unstructured_limiter_min_pressure_filter     = true

-- viscous settings
config.viscous = true
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"

-- Set temporal integration settings
NewtonKrylovGlobalConfig{
   -- phases
   number_of_phases = 3,
   max_steps_in_initial_phases = { 500, 500 },
   phase_changes_at_relative_residual = { 1.0e-03, 1.0e-06 },

   -- Preconditioner
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,

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

   -- output setting
   number_of_steps_for_setting_reference_residuals = 0,
   steps_between_status = 1,
   write_loads = true,
   steps_between_loads_update = 1,
   total_snapshots = 5,
   steps_between_snapshots = 1,
   steps_between_diagnostics = 1,
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
   start_cfl = 1.0,
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
