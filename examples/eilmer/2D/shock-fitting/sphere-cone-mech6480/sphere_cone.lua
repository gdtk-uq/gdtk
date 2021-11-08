-- sphere_cone.lua
-- Sphere Cone in ideal air flow with shock-fitting boundary
-- MECH6480 class example.
-- PJ 2021-10-21
config.title = "Sphere-cone in ideal air flow."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
setGasModel('ideal-air-gas-model.lua')
p_inflow = 1.0e3 -- Pa
T_inflow = 200.0 -- K
v_inflow = 7.0*math.sqrt(1.4*287.1*T_inflow) -- m/s
print("inflow:", p_inflow, "Pa,", T_inflow, "K,", v_inflow, "m/s")
initial = FlowState:new{p=p_inflow/10.0, T=300.0}
inflow = FlowState:new{p=p_inflow, T=T_inflow, velx=v_inflow}
--
print("Building grid.")
R = 0.010 -- m
L = 0.030 -- m
theta = 10.0*math.pi/180.0
sn = math.sin(theta); cs = math.cos(theta)
--
a1 = {x=0.0, y=0.0}; a2 = {x=-R, y=0.0}
b1 = {x=R-R*sn, y=R*cs}; b2 = {x=R-2*R*sn, y=2*R*cs}
c = {x=R, y=0.0}
d1 = {x=L, y=b1.y+(L-R*cs)*sn/cs}; d2 = {x=L, y=d1.y+1.5*R}
--
a1b1 = Arc:new{p0=a1, p1=b1, centre=c}
b1d1 = Line:new{p0=b1, p1=d1}
a2b2 = Arc:new{p0=a2, p1=b2, centre=c}
b2d2 = Line:new{p0=b2, p1=d2}
a2a1 = Line:new{p0=a2, p1=a1}
b2b1 = Line:new{p0=b2, p1=b1}
d2d1 = Line:new{p0=d2, p1=d1}
--
psurf0 = makePatch{north=b2b1, east=a1b1, south=a2a1, west=a2b2}
grid0 = StructuredGrid:new{psurface=psurf0, niv=20, njv=40}
psurf1 = makePatch{north=d2d1, east=b1d1, south=b2b1, west=b2d2}
grid1 = StructuredGrid:new{psurface=psurf1, niv=20, njv=60}
-- To allow the fitted shock to be the west edge of a single FluidBlock Array
-- we need to start with a single grid that will be subdivided.
grid0:joinGrid(grid1, "north")
--
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
