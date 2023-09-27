-- corner.lua
--
-- Mach 2.0 flow expanding around a convex corner,
-- details are taken from the AIAA verification project.
--
-- Kyle A. Damm [2022-08-01]
--

-- Try to get the file format from environment variable
fileFmt = os.getenv("LMR_FILE_FORMAT") or "rawbinary"

job_title = "Mach 2 air flowing around a convex corner."
config.title = job_title
print(job_title)
print("File format in use: ", fileFmt)

-- General settings
config.dimensions = 2
config.axisymmetric = false
config.print_count = 1
config.new_flow_format = true
config.save_residual_values = false
config.save_limiter_values = false
config.flow_format = fileFmt
config.grid_format = fileFmt

-- ==========================================================
-- Define the flow domain using a native grid.
-- ==========================================================

-- geometric parameters
L = 1.0
theta = 10.0 * (math.pi/180.0)
dy = L*math.tan(theta)

-- define vertices
p0 = Vector3:new{x=-L,  y=0.0}
p1 = Vector3:new{x=0.0, y=0.0}
p2 = Vector3:new{x=L,   y=-dy}
p3 = Vector3:new{x=L,   y=L-dy}
p4 = Vector3:new{x=0.0, y=L}
p5 = Vector3:new{x=-L,  y=L}

-- define paths
l01 = Line:new{p0=p0, p1=p1}
l12 = Line:new{p0=p1, p1=p2}
l02 = Polyline:new{segments={l01,l12}}
l23 = Line:new{p0=p2, p1=p3}
l54 = Line:new{p0=p5, p1=p4}
l43 = Line:new{p0=p4, p1=p3}
l53 = Polyline:new{segments={l54,l43}}
l43 = Line:new{p0=p4, p1=p3}
l05 = Line:new{p0=p0, p1=p5}

-- define patch
quad0 = makePatch{north=l53, east=l23, south=l02, west=l05}

-- define grid
nx = 30
ny = 15
grid0 = registerGrid{
   grid=StructuredGrid:new{psurface=quad0, niv=nx+1, njv=ny+1},
   fsTag="initial",
   bcTags={north="outflow",east="outflow",west="inflow"}
}

-- ==========================================================
-- Freestream conditions
-- ==========================================================
nsp, nmodes, gm = setGasModel('ideal-air.gas')
gs = GasState:new{gm}
M_inf = 2.0
gs.p  = 100.0e3 -- Pa
gs.T  = 300.0 -- K
gm:updateThermoFromPT(gs)
gm:updateSoundSpeed(gs)
V_inf = M_inf*gs.a
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}
initial = inflow

-- ==========================================================
-- Block definitions
-- ==========================================================
flowDict = {}
flowDict["initial"] = initial
flowDict["inflow"] = inflow

bcDict = {
   inflow = InFlowBC_Supersonic:new{flowState=inflow},
   outflow = OutFlowBC_Simple:new{}
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
config.apply_limiter = false

-- viscous flux settings
config.viscous = false

-- Set temporal integration settings
config.residual_smoothing = false

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 50,
   stop_on_relative_residual = 1.0e-12,
   number_of_phases = 1,
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = false,
   max_linear_solver_iterations = 50,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
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
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 5,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.99,
   start_cfl = 1.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}


