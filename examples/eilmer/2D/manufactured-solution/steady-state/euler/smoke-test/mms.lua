-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid
--          2015-06-03
--          Added scripting to configure case based on 'case.txt'
--
-- Ported to steady-state solver (Kyle A. Damm, 2018)

config.title = "Method of Manufactured Solutions."
print(config.title)
config.dimensions = 2

file = io.open("case.txt", "r")
scale = tonumber(file:read("*line"))
fluxCalc = tostring(file:read("*line"))
derivCalc = tostring(file:read("*line"))
xOrder = tonumber(file:read("*line"))
blocking = tostring(file:read("*line"))
ncells = tonumber(file:read("*line"))
file:close()

setGasModel('very-viscous-air.lua')
R = 287.0
p0 = 1.0e5; T0 = p0/R;
u0 = 800.0; v0 = 800.0
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}

p00 = Vector3:new{x=0.0, y=0.0}
p10 = Vector3:new{x=1.0, y=0.0}
p01 = Vector3:new{x=0.0, y=1.0}
p11 = Vector3:new{x=1.0, y=1.0}
nicell = ncells; njcell = ncells
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
			  niv=nicell+1, njv=njcell+1}

bcList = {}
bcList['north'] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList['east'] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList['south'] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList['west'] = UserDefinedBC:new{fileName='udf-bc.lua'}
config.apply_bcs_in_parallel = false

if blocking == 'single' then
   blks = FluidBlock:new{grid=grid, fillCondition=initial, bcList=bcList, label="blk"}
   SBlock2UBlock(fluidBlocks[1])
else
   blks = FBArray:new{grid=grid, fillCondition=initial, bcList=bcList,
		      nib=4, njb=4, label="blk"}
   for i=1,16 do
      SBlock2UBlock(fluidBlocks[i])
   end
end

config.include_ghost_cells_in_spatial_deriv_clouds = true
config.interpolation_order = xOrder
config.flux_calculator = fluxCalc
config.spatial_deriv_calc = derivCalc
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.viscous = false
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false
-- Set simulation parameters
config.dt_loads = 1.0e-3
config.print_count = 100
SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 100;
   start_preconditioning = 1,

   use_scaling = true,
   use_complex_matvec_eval = true,

   number_pre_steps = 10,
   number_total_steps = 100000,
   stop_on_relative_global_residual = 1.0e-12,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 10,
   max_restarts = 10,

   -- Settings for start-up phase
   number_start_up_steps = 100,
   cfl0 = 1.0,
   eta0 = 0.1,
   tau0 = 1.0,
   sigma0 = 1.0e-30,
   p0 = 1.0,

   -- Settings for inexact Newton phase
   cfl1 = 100.0,
   tau1 = 1.0,
   sigma1 = 1.0e-30,
   p1 = 1.0,
   eta1 = 0.1,
   eta1_min = 0.1,
   eta_strategy = "geometric",

   -- Settings control write-out
   snapshots_count = 50,
   number_total_snapshots = 5,
   write_diagnostics_count = 20
}
