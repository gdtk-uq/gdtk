-- cyl-sf.lua
-- Shock fitting boundary condition example simulation:
-- Cylinder in ideal air flow
-- Kyle Damm Jan 2015

job_title = "Cylinder in ideal air flow with shock fitting boundary."
print(job_title)

config.dimensions = 2
config.title = job_title

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=100.0e3/3.0, T=200.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=100.0e3, T=300.0, velx=2430.0, vely=0.0}
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)

print "Building grid."
a = Vector3:new{x=-3.0, y=0.0}
b = Vector3:new{x=-1.0, y=0.0}
c = Vector3:new{ x=0.0, y=3.0}
d = Vector3:new{ x=0.0, y=1.0}
p = Vector3:new{ x=0.0, y=0.0}

ac = Arc:new{p0=a, p1=c, centre=p}
bd = Arc:new{p0=b, p1=d, centre=p}
ab = Line:new{p0=a, p1=b}
cd = Line:new{p0=c, p1=d}

psurf = makePatch{north=cd, east=bd, south=ab, west=ac}
grid = StructuredGrid:new{psurface=psurf, niv=25, njv=30}

-- We can leave east and south as slip-walls
blk0 = SBlock:new{grid=grid, fillCondition=initial}
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_ShockFitting:new{flowCondition=inflow, label="shockfitting-inflow-boundary"}
blk0.bcList[north] =  OutFlowBC_Simple:new{label="outflow-boundary"}

-- Set a few more config options
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "euler"
config.max_time = 15.0e-3
config.max_step = 100000
config.dt_init = 1.0e-9
config.cfl_value = 0.25 
config.dt_plot = 5.0e-4

-- moving grid flag
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 1.383e-03
