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
config.gasdynamic_update_scheme = 'rkl2'
config.max_time                 = 1.75*0.001 - 1.0e-08 -- to allow for round-off when trying to hit 0.001 
config.dt_init                  = dt
config.fixed_time_step          = true
config.max_step                 = 1e6
config.max_attempts_for_step    = 1
config.dt_plot                  = config.max_time/20

