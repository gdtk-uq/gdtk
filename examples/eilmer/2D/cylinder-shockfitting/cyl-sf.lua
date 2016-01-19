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
a = Vector3:new{-3.0,    0.0}
b = Vector3:new{-1.0,    0.0}
c = Vector3:new{0.0,     3.0}
d = Vector3:new{0.0,     1.0}
p = Vector3:new{0.0,     0.0}

ac = Arc:new{a, c, p}
bd = Arc:new{b, d, p}
ab = Line:new{a, b}
cd = Line:new{c, d}

psurf = makePatch{cd, bd, ab, ac}
grid = StructuredGrid:new{psurface=psurf, niv=25, njv=30}

-- We can leave east and south as slip-walls
blk0 = SBlock:new{grid=grid, fillCondition=initial}
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_ShockFitting:new{flowCondition=inflow, label="shockfitting-inflow-boundary"}
blk0.bcList[north] =  OutFlowBC_Simple:new{label="outflow-boundary"}

-- Set a few more config options
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "euler"
config.max_time = 5.0e-2
config.max_step = 400000
config.dt_init = 1.0e-9
config.cfl_value = 0.25 
config.dt_plot = 5.0e-4

-- moving grid flag
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 1.383e-03
