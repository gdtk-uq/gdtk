--
-- Method of Manufactored Solutions, temporal order of accuracy verification example
-- 
-- author: Kyle A. Damm
-- date:   29-08-2023
--

job_title = "Method of Manufactured Solutions, temporal verification case."
config.title = job_title
print(job_title)

-- ==========================================================
-- General settings
-- ==========================================================
config.print_count = 1
-- Parse dt value from file
file = io.open("case.txt", "r")
dt = tonumber(file:read("*line"))
file:close()
-- read constants from file
dofile('constants.txt')

-- ==========================================================
-- Flow conditions
-- ==========================================================
setGasModel('very-viscous-air.lua')
-- Load initial condition fill functions from file
dofile('fill-fn.lua')

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
grid = StructuredGrid:new{psurface=patch, niv=nx0+1, njv=ny0+1}

-- define structured fluidblock
blk = FluidBlock:new{grid=grid, initialState=gasFillFn, label="blk0"}

-- define boundary conditions
config.apply_bcs_in_parallel = false
blk.bcList[north] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                         UpdateThermoTransCoeffs:new()
   }
}
blk.bcList[east] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                         UpdateThermoTransCoeffs:new()
   }
}
blk.bcList[south] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                         UpdateThermoTransCoeffs:new()
   }
}
blk.bcList[west] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                         UpdateThermoTransCoeffs:new()
   }
}

identifyBlockConnections()

-- ==========================================================
-- Solver settings
-- ==========================================================

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
config.max_time                 = 1.75*0.001 - 1.0e-08 -- to allow for round-off when trying to hit 0.001 
config.dt_init                  = dt
config.dt_plot                  = config.max_time/20
config.fixed_time_step          = true
config.with_local_time_stepping = true
SteadyStateSolver{
   -- transient or steady-state
   temporal_integration_mode = 2,

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
   physicality_check_theta = 0.5,
   use_complex_matvec_eval = true,

   -- freeze limiter
   limiter_freezing_residual_reduction = 1.0e-20,
   limiter_freezing_count = 1,

   -- stopping criteria
   number_pre_steps = 0,
   number_total_steps = 5e6,
   stop_on_relative_global_residual = 1.0e-12,
   stop_on_absolute_global_residual = 1.0e-03,
   
   -- settings for GMRES iterative solver
   max_outer_iterations = 100,
   max_restarts = 0,

   -- CFL settings
   residual_based_cfl_scheduling = true,
   inviscid_cfl = true,

   -- start-up phase settings
   number_start_up_steps = 0,
   sigma0 = 1.0e-50,
   cfl0 = 1.0,

   -- final phase settings
   sigma1 = 1.0e-50,
   cfl1 = 1.0e6,
   p1 = 1.0,
   tau1 = 1.0e-16,
   eta1 = 1.0e-04,
   eta_strategy = "constant",

   -- output settings
   snapshots_count = 1,
   number_total_snapshots = 5000,
   write_diagnostics_count = 1,
   write_loads_count = 250,
}
