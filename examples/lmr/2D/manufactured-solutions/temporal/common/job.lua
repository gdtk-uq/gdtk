--
-- Method of Manufactored Solutions, temporal order of accuracy verification example
-- 
-- author: Kyle A. Damm
-- date:   29-08-2023
-- ported to LMR by NNG, May 2025

job_title = "Method of Manufactured Solutions, temporal verification case."
print(job_title)

-- ==========================================================
-- General settings
-- ==========================================================
config.print_count = 1
-- Parse dt and gas_dynamic value from file
dofile('config.txt')
-- read constants from file
dofile('constants.txt')

-- ==========================================================
-- Flow conditions
-- ==========================================================
setGasModel('very-viscous-air.lua')
-- Load initial condition fill functions from file
dofile('fill-fn.lua')
flowStates = {initial=gasFillFn}

-- ==========================================================
-- Define the flow domain using an native grid
-- ==========================================================
config.dimensions = 2
ncells = 10
nx0 = ncells
ny0 = ncells

-- define nodes
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=L, y=0.0}
c = Vector3:new{x=0.0, y=L}
d = Vector3:new{x=L, y=L}

-- define patch
patch = makePatch{north=Line:new{p0=c, p1=d},
                  east=Line:new{p0=b, p1=d},
                  south=Line:new{p0=a, p1=b},
                  west=Line:new{p0=a, p1=c}}

-- define grid of cells
sgrid = StructuredGrid:new{psurface=patch, niv=nx0+1, njv=ny0+1}
grid = registerFluidGrid{
   grid=UnstructuredGrid:new{sgrid=sgrid},
   fsTag="initial",
   bcTags={[Face.north]="udf_north", [Face.east]="udf_east", [Face.south]="udf_south", [Face.west]="udf_west"}

}


-- define boundary conditions
bcs = {
   udf_north = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                            UpdateThermoTransCoeffs:new()
      }
   },
   udf_east = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                            UpdateThermoTransCoeffs:new()
      }
   },
   udf_south = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                            UpdateThermoTransCoeffs:new()
      }
   },
   udf_west = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                            UpdateThermoTransCoeffs:new()
      }
   }
}

makeFluidBlocks(bcs, flowStates)

-- ==========================================================
-- Solver settings
-- ==========================================================
config.solver_mode = solver_mode
config.apply_bcs_in_parallel = false

-- invsicid flux settings
config.flux_calculator     = 'ldfss2'
config.interpolation_order = 2
config.apply_limiter       = false
config.extrema_clipping    = false

-- viscous flux settings
config.viscous            = true
config.spatial_deriv_calc = 'least_squares'

-- source term settings
config.udf_source_terms                    = true
config.udf_source_terms_file               = 'udf-source-terms.lua'
config.eval_udf_source_terms_at_each_stage = true

-- Set temporal integration settings
if config.solver_mode == 'dual_time_stepping' then
   config.dualtimestepping_update_scheme = update_scheme
else
   config.gasdynamic_update_scheme = update_scheme
end
config.max_time                 = 1.75*0.001 - 1.0e-08 -- to allow for round-off when trying to hit 0.001 
config.dt_init                  = dt
config.fixed_time_step          = true
config.max_step                 = 1e6
config.max_attempts_for_step    = 1
config.dt_plot                  = config.max_time/20

-- These are settings from Rowan's LMR Navier-stokes MMS
config.include_ghost_cells_in_spatial_deriv_clouds = true

-- These are settings used only for dual time stepping schemes
NewtonKrylovGlobalConfig{
   -- phases
   number_of_phases = 1,

   -- preconditioner settings
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,

   -- nonlinear system solver settings
   max_newton_steps = 100,
   max_consecutive_bad_steps = 10,
   stop_on_relative_residual = 1.0e-10,

   -- linear system solver settings
   frechet_derivative_perturbation = 1.0e-50,
   use_scaling = true,
   max_linear_solver_iterations = 100,
   max_linear_solver_restarts = 0,

   -- continuation settings
   inviscid_cfl_only = true,
   use_line_search = false,
   use_residual_smoothing = false,
   use_physicality_check = true,
   allowable_relative_mass_change = 0.9,
   min_relaxation_factor_for_update = 0.1,
   min_relaxation_factor_for_cfl_growth = 0.5,
}

NewtonKrylovPhase:new{
   -- preconditioner settings
   frozen_preconditioner = true,
   use_adaptive_preconditioner = true,
   steps_between_preconditioner_update = 5,

   -- linear system solver settings
   linear_solve_tolerance = 0.001,
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 1,

   -- cfl settings
   use_auto_cfl = true,
   use_local_timestep = true,
   threshold_relative_residual_for_cfl_growth = 0.99,
   start_cfl = 1.0e+01,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}

