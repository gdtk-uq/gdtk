-- An Eilmer job script to simulate a Busemman diffuser.
--
-- Example taken from Moelder (2003)
--
-- Reference:
-- 
-- Moelder, S. (2003)
-- A Benchmark for Internal Flow CFD Codes
-- Computational Fluid Dynamics Journal, 12(2):47, pp.408--414
--
-- Related work in collection:
-- There is a version of this simulation built for Puffin.
-- It is located at: gdtk/examples/puffin/busemann-diffuser
--
-- Author: RJG
-- Date: 2024-07-06

contourFile = "bd-contour.dat"
-- Load the diffuser properties into working space of script
dofile("bd-props.txt")
print(string.format("Busemann diffuser with M1=%.2f, M2=%.2f, M3=%.2f", M1, M2, M3))

config.dimensions = 2
config.axisymmetric = true
config.solver_mode = "steady"

-- Set up geometry and grid.
--
--   ._________          
-- e |         \_______  bodyside  
-- n |                 \_____
-- t |                       \____
-- r |                            \__      thrt
-- a |   gridArray[0]                \.______________. e
-- n |                                |              | x
-- c |                                | gridArray[1] | i
-- e |                                |              | t
--   .___ _ ___ _ ___ _ ___ _ ___ _ __.__ ___ _ ___ _._ ___ _ CL
--            symmetry0                   symmetry1
--  
thrtLength = 0.2
spline = Spline2:new{filename=contourFile}
bodyside = ArcLengthParameterizedPath:new{underlying_path=spline}
thrtStart = bodyside(1.0)
thrt = Line:new{p0={x=thrtStart.x, y=thrtStart.y},
                p1={x=thrtStart.x + thrtLength, y=thrtStart.y}}
symmetry0 = Line:new{p0={x=bodyside(0.0).x, y=0.0}, p1={x=thrtStart.x, y=0.0}}
symmetry1 = Line:new{p0={x=thrtStart.x, y=0.0}, p1={x=thrt(1.0).x, y=0.0}}
entrance = Line:new{p0=symmetry0(0.0), p1=bodyside(0.0)}
thrtPlane = Line:new{p0=symmetry0(1.0), p1=bodyside(1.0)}
exit = Line:new{p0=symmetry1(1.0), p1=thrt(1.0)}

-- patches
quad0 = ControlPointPatch:new{north=bodyside, south=symmetry0, east=thrtPlane, west=entrance,
                              ncpi=10, ncpj=4, guide_patch="channel"}
--quad0 = CoonsPatch:new{north=bodyside, south=symmetry0, east=thrtPlane, west=entrance}
quad1 = CoonsPatch:new{north=thrt, south=symmetry1, east=exit, west=thrtPlane}
-- gridarrays
-- Let's go for squarish cells at entrance, and then expect them to be squashed in
-- aspect ratio as we head towards the throat section.
nycells = 20
nxcells0 = math.ceil(nycells*symmetry0:length()/entrance:length())
nxcells1 = math.ceil(nxcells0*symmetry1:length()/symmetry0:length())
print("nxcells0= ", nxcells0, " nxcells1= ", nxcells1)
njb = 1
nib1 = 1
nib0 = math.ceil(nib1*symmetry0:length()/symmetry1:length())
print("nib0= ", nib0)
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=quad0, niv=nxcells0+1, njv=nycells+1},
   nib=nib0, njb=njb,
   fsTag='initial',
   bcTags={west='inflow'}
}
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=quad1, niv=nxcells1+1, njv=nycells+1},
   nib=nib1, njb=njb,
   fsTag='initial',
   bcTags={east='outflow'}
}
identifyGridConnections()

-- Set up some flow conditions
_, __, gm = setGasModel("ideal-air.lua")
gs = GasState:new{gm}
q1 = 50e3 -- Pa
T1 = 250.0 -- K
gs.T = T1
gamma = gm:gamma(gs)
gs.p = 2*q1/((gamma-1)*M1*M1)
gm:updateThermoFromPT(gs)
gm:updateSoundSpeed(gs)
V1 = M1 * gs.a
print("V1= ", V1)
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V1, vely=0.0}

flowDict = {
   initial=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
makeFluidBlocks(bcDict, flowDict)

mpiDistributeBlocks{ntasks=4}

-- Simulation parameters
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.extrema_clipping = false

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-10,
   max_newton_steps = 500,
   number_of_phases = 1,
   max_linear_solver_iterations = 50,
   total_snapshots = 3,
   steps_between_status = 1,
   steps_between_snapshots = 10,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   use_auto_cfl = true,
   start_cfl = 1.0,
   threshold_relative_residual_for_cfl_growth = 0.9,
   auto_cfl_exponent = 0.9
}














