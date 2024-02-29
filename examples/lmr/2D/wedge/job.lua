-- Mach 3.0 flow over a 15 degree wedge,
-- details are taken from the AIAA verification project.
--
-- Reference:
-- Ghia et al. (2010)
-- The AIAA code verification project-test cases for CFD code verification,
-- AIAA Paper 2010-0125
--
-- Kyle A. Damm [2022-08-01]
--
-- Updated for Eilmer v5 by RJG, 2023-07-05
--

-- Try to get the file format from environment variable
fileFmt = os.getenv("LMR_FILE_FORMAT") or "rawbinary"

job_title = "Mach 3 air flowing over a 15 degree wedge."
print(job_title)

-- General settings
config.solver_mode = 'steady'
config.dimensions = 2
config.axisymmetric = false
config.print_count = 1
config.save_residual_values = true
config.save_limiter_values = true
config.field_format = fileFmt
config.grid_format = fileFmt

-- ==========================================================
-- Freestream conditions
-- ==========================================================
nsp, nmodes, gm = setGasModel('ideal-air.gas')
gs = GasState:new{gm}
M_inf = 3.0
gs.p  = 100.0e3 -- Pa
gs.T  = 300.0 -- K
gm:updateThermoFromPT(gs)
gm:updateSoundSpeed(gs)
V_inf = M_inf*gs.a
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}
initial = inflow

-- ==========================================================
-- Define the flow domain using a native grid.
-- ==========================================================

-- geometric parameters
L = 1.0
theta = 15.0 * (math.pi/180.0)
dy = L*math.tan(theta)

-- define vertices
p0 = {x=-L,  y=0.0}
p1 = {x=0.0, y=0.0}
p2 = {x=L,   y=dy}
p3 = {x=L,   y=L}
p4 = {x=0.0, y=L} 
p5 = {x=-L,  y=L}

-- define paths
l01 = Line:new{p0=p0, p1=p1}
l12 = Line:new{p0=p1, p1=p2}
l23 = Line:new{p0=p2, p1=p3}
l43 = Line:new{p0=p4, p1=p3}
l54 = Line:new{p0=p5, p1=p4}
l05 = Line:new{p0=p0, p1=p5}
l14 = Line:new{p0=p1, p1=p4}

-- define patch
quad0 = makePatch{north=l54, east=l14, south=l01, west=l05}
quad1 = makePatch{north=l43, east=l23, south=l12, west=l14}

-- define grid
nx0 = 15
nx1 = nx0
ny = 15
-- Build a structured grids, join them and then 
-- convert to unstructured when we register the grid.
sgrid0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1}
sgrid1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1}
sgrid0:joinGrid(sgrid1, "east")

grid0 = registerFluidGrid{
   grid=UnstructuredGrid:new{sgrid=sgrid0},
   fsTag="initial",
   bcTags={[Face.west]="inflow", [Face.east]="outflow", [Face.south]="wall", [Face.north]="outflow"}
}

flowDict = {}
flowDict["initial"] = initial

bcDict = {
   inflow = InFlowBC_Supersonic:new{flowState=inflow},
   outflow = OutFlowBC_SimpleExtrapolate:new{},
   wall = WallBC_WithSlip:new{}
}

makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver configuration
-- ==========================================================

-- invsicid flux settings
config.flux_calculator= "ausmdv"
config.apply_entropy_fix = false
config.interpolation_order = 2
config.extrema_clipping = false
config.thermo_interpolator = "rhop"
config.unstructured_limiter = "park"

-- viscous flux settings
config.viscous = false

-- Set temporal integration settings
config.residual_smoothing = false

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 1500,
   stop_on_relative_residual = 1.0e-14,
   number_of_phases = 1,
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = false,
   use_scaling = true,
   max_linear_solver_iterations = 50,
   max_linear_solver_restarts = 0,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 3,
   steps_between_snapshots = 5,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   frozen_preconditioner = true,
   steps_between_preconditioner_update = 5,
   use_adaptive_preconditioner = false,
   frozen_limiter_for_jacobian = false,
   linear_solve_tolerance = 0.01,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 1.0,
   start_cfl = 1.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}
