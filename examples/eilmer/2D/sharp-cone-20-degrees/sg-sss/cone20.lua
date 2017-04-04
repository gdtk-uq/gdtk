-- cone20.lua
-- Simple job-specification file for e4prep -- for use with Eilmer4
-- PJ & RG
-- 2015-02-24 -- adapted from the Python version of cone20
-- 2017-04-05 -- this version uses the eilmer4 steady-state solver.

-- We can set individual attributes of the global data object.
config.title = "Mach 1.5 flow over a 20 degree cone."
print(config.title)
config.dimensions = 2
config.axisymmetric = true

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}

-- Demo: Verify Mach number of inflow and compute dynamic pressure.
Q = inflow:toTable()
print("T=", Q.T, "density=", Q.rho, "sound speed= ", Q.a)
print("inflow Mach number=", 1000.0/Q.a)
print("dynamic pressure q=", 1/2*Q.rho*1.0e6)

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.2, y=0.0}
c = Vector3:new{x=1.0, y=0.29118}
d = Vector3:new{x=1.0, y=1.0}
e = Vector3:new{x=0.2, y=1.0}
f = Vector3:new{x=0.0, y=1.0}
ab = Line:new{p0=a, p1=b} -- lower boundary, axis
bc = Line:new{p0=b, p1=c} -- lower boundary, cone surface
fe = Line:new{p0=f, p1=e}; ed = Line:new{p0=e, p1=d} -- upper boundary
af = Line:new{p0=a, p1=f} -- vertical line, inflow
be = Line:new{p0=b, p1=e} -- vertical line, between quads
cd = Line:new{p0=c, p1=d} -- vertical line, outflow
quad0 = makePatch{north=fe, east=be, south=ab, west=af}
quad1 = makePatch{north=ed, east=cd, south=bc, west=be, gridType="ao"}
-- Mesh the patches, with particular discretisation.
nx0 = 10; nx1 = 30; ny = 40
grid0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1}
-- Define the flow-solution blocks.
blk0 = SBlock:new{grid=grid0, fillCondition=inflow}
blk1 = SBlock:new{grid=grid1, fillCondition=inflow}
-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_Supersonic:new{flowCondition=inflow}
blk1.bcList[east] = OutFlowBC_Simple:new{}

-- add history point 1/3 along length of cone surface
setHistoryPoint{x=2*b.x/3+c.x/3, y=2*b.y/3+c.y/3}
-- add history point 2/3 along length of cone surface
setHistoryPoint{ib=1, i=math.floor(2*nx1/3), j=0}

-- Do a little more setting of global data.
config.print_count = 1

SteadyStateSolver{
   use_preconditioning = false,
   number_pre_steps = 3,
   number_total_steps = 120,
   -- Settings for FGMRES iterative solver
   max_outer_iterations = 10,
   max_restarts = 0,
   -- Settings for start-up phase
   number_start_up_steps = 10,
   cfl0 = 2.0,
   eta0 = 0.5,
   tau0 = 0.1,
   sigma0 = 5.0e-6,
   -- Settings for inexact Newton phase
   cfl1 = 100000.0,
   tau1 = 0.1,
   sigma1 = 5.0e-6,
   eta1 = 0.0001,
   eta_strategy = "constant",
   -- Settings control write-out
   snapshots_count = 20,
   number_total_snapshots = 1,
   write_diagnostics_count = 1
}

