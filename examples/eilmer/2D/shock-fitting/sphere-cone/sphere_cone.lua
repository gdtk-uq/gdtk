-- sphere_cone.lua
-- Sphere Cone in ideal air flow with shock-fitting boundary
-- KD  2016-12-16
-- PJ  2021-08-07 Rebuild to use a single FluidBlock Array.
--
config.title = "Sphere-cone in ideal air flow."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
setGasModel('ideal-air-gas-model.lua')
p_inflow = 276.746 -- Pa
T_inflow = 46.2227 -- K
v_inflow = 1068.436 -- m/s
initial = FlowState:new{p=p_inflow/3.0, T=T_inflow, velx=0.0, vely=0.0}
inflow = FlowState:new{p=p_inflow, T=T_inflow, velx=v_inflow, vely=0.0}
--
print "Building grid."
R = 0.00508 -- m
L = 0.18476 -- m
RL = 0.03683 -- m
--
p = {x=0.0, y=0.0}
a = {x=-1*R, y=0.0}
b = {x=-R*math.cos(80.0 * math.pi/180.0), y=R*math.sin(80.0 * math.pi/180.0)}
c = {x=L, y=RL}
e = {x=-3*R, y=0}
f = {x=-3*R*math.cos(80.0 * math.pi/180.0), y=3*R*math.sin(80.0 * math.pi/180.0)}
g = {x=L, y=RL+8*R}
--
ab = Arc:new{p0=a, p1=b, centre=p}
bc = Line:new{p0=b, p1=c}
ef = Arc:new{p0=e, p1=f, centre=p}
fg = Line:new{p0=f, p1=g}
ea = Line:new{p0=e, p1=a}
fb = Line:new{p0=f, p1=b}
gc = Line:new{p0=g, p1=c}
--
-- cluster_y = RobertsFunction:new{end0=true, end1=false, beta=1.05}
cluster_x = RobertsFunction:new{end0=true, end1=false, beta=1.04}
cluster_none = RobertsFunction:new{end0=false, end1=false, beta=1.05}
cflist0 = {north=cluster_none, east=cluster_none, south=cluster_none, west=cluster_none}
cflist1 = {north=cluster_none, east=cluster_x, south=cluster_none, west=cluster_x}
--
psurf0 = makePatch{north=fb, east=ab, south=ea, west=ef}
grid0 = StructuredGrid:new{psurface=psurf0, cfList=cflist0, niv=60, njv=50}
psurf1 = makePatch{north=gc, east=bc, south=fb, west=fg}
grid1 = StructuredGrid:new{psurface=psurf1, cfList=cflist1, niv=60, njv=100}
--
-- To allow the fitted shock to be the west edge of a single FluidBlock Array
-- we need to start with a single grid that will be subdivided.
grid0:joinGrid(grid1, "north")
blk = FBArray:new{grid=grid0, initialState=initial,
                  bcList={west=InFlowBC_ShockFitting:new{flowState=inflow},
                          north=OutFlowBC_Simple:new{}},
                  nib=1, njb=4}
identifyBlockConnections()
--
-- Set a few more config options
body_flow_time = L/v_inflow
print("body_flow_time=", body_flow_time)
config.gasdynamic_update_scheme = "backward_euler"
config.max_time = 10*body_flow_time
config.max_step = 400000
config.dt_init = 1.0e-8
config.cfl_value = 0.5
config.dt_plot = body_flow_time
-- We want the starting pulse to clear the flow domain before allowing
-- the inflow boundary to move down into the shock location.
-- Until then, the shock is too complicated.
config.shock_fitting_delay = 2*body_flow_time
config.grid_motion = "shock_fitting"
config.max_invalid_cells = 10
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
