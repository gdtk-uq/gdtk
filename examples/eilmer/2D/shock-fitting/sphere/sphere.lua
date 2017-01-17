-- sphere.lua -- Sphere in ideal air flow with shock fitting boundary
-- KD  2016-01-01

job_title = "Sphere in ideal air flow."
print(job_title)

config.dimensions = 2
config.title = job_title
config.axisymmetric = true
setGasModel('ideal-air-gas-model.lua')
p_inflow = 276.746 -- Pa
T_inflow = 46.2227 -- K
v_inflow = 1068.436 -- m/s
initial = FlowState:new{p=p_inflow/3.0, T=T_inflow, velx=0.0, vely=0.0}
inflow = FlowState:new{p=p_inflow, T=T_inflow, velx=v_inflow, vely=0.0}
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
print "Building grid."
R = 0.0508 -- m

a = Vector3:new{x=-3*R,    y=0.0}
b = Vector3:new{x=-1*R,    y=0.0}
c = Vector3:new{x=0.0,     y=3*R}
d = Vector3:new{x=0.0,     y=1*R}
p = Vector3:new{x=0.0,     y=0.0}

ac = Arc:new{p0=a, p1=c, centre=p}
bd = Arc:new{p0=b, p1=d, centre=p}
ab = Line:new{p0=a, p1=b}
cd = Line:new{p0=c, p1=d}

psurf = makePatch{north=cd, east=bd, south=ab, west=ac}
grid = StructuredGrid:new{psurface=psurf, niv=100, njv=100}

-- We can leave east and south as slip-walls
blk = SBlockArray{grid=grid, fillCondition=initial,
		  bcList={west=InFlowBC_ShockFitting:new{flowCondition=inflow},
			  north=OutFlowBC_Simple:new{}},
		  nib=1, njb=8}
identifyBlockConnections()

-- Set a few more config options
body_flow_time = 2*R/v_inflow
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.max_step = 400000
config.dt_init = 1.0e-9
config.cfl_value = 0.5 
config.dt_plot = body_flow_time
config.interpolation_order = 2
-- moving grid flag
body_flow_time = 2*R/v_inflow
config.shock_fitting_delay = 1*body_flow_time
config.grid_motion = "shock_fitting"
config.max_time = 25*body_flow_time
