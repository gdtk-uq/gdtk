--MHD Shock tube example
--Classic MHD validation case is the Brio and Wu test case
--Alternate test case is the Ryu and Jones test case

--set attributes of the simulation
config.title = "MHDShockTube"
print(config.title)
config.dimensions = 2
config.flux_calculator = "hlle"
config.MHD = true
config.divergence_cleaning = false
config.axisymmetric = false
config.stringent_cfl = true

--set gas model to ideal air
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')

--Define the left and right states
R = 287

--Sod's Shock Tube- classic validation case for any CFD code
p_left = 1.0e4
p_right = 0.1e4
rho_left = 1
rho_right = 0.125
Bx = 0
By_left = 0
By_right = 0
T_left = p_left / (rho_left * R)
T_right = p_right / (rho_right * R)
left = FlowState:new{p=p_left, T=T_left, Bx=Bx, By=By_left}
right = FlowState:new{p=p_right, T=T_right, Bx=Bx, By=By_right}

--[[Brio-Wu Shock Tube- basic MHD validation case
p_left = 1.0 / (4e-4 * 4e-4 * 8e-2)
p_right = 0.1 / (4e-4 * 4e-4 * 8e-2)
rho_left = 1 / (8e-2 * 8e-2 * 8e-2 )
rho_right = 0.125 / (8e-2 * 8e-2 * 8e-2)
Bx = 0.75 / (4e-4 * math.sqrt(8e-2))
By_left = 1.0 / (4e-4 * math.sqrt(8e-2))
By_right = -1.0 / (4e-4 * math.sqrt(8e-2))
T_left = p_left / (rho_left * R)
T_right = p_right / (rho_right * R)
left = FlowState:new{p=p_left, T=T_left, Bx=Bx, By=By_left}
right = FlowState:new{p=p_right, T=T_right, Bx=Bx, By=By_right}]]

--[[For the Ryu and Jones Test 2A
Bx = 2e2/math.sqrt(4*math.pi)
rho_left = 1.08
vx_left = 1.2e2
vy_left = 0.01e2
vz_left = 0.5e2
By_left = 3.6e2/math.sqrt(4*math.pi)
Bz_left = 2e2/math.sqrt(4*math.pi)
p_left = 0.95e4
T_left = p_left / (R * rho_left)
left = FlowState:new{p=p_left, T=T_left, velx=vx_left, vely=vy_left, velz=vz_left, Bx=Bx, By=By_left, Bz=Bz_left}

rho_right = 1
vx_right = 0e2
vy_right = 0e2
vz_right = 0e2
By_right = 4e2/math.sqrt(4*math.pi)
Bz_right = 2e2/math.sqrt(4*math.pi)
p_right = 1e4
T_right = p_right / (R * rho_right)
right = FlowState:new{p=p_right, T=T_right, velx=vx_right, vely=vy_right, velz=vz_right, Bx=Bx, By=By_right, Bz=Bz_right}]]

--Define the simulation domain; set to a unit length
x_start = 0.0
x_mid = 0.5
x_end = 1.0
y_bot = 0.0
y_top = 0.5

a = Vector3:new{x=x_start, y=y_bot}
b = Vector3:new{x=x_mid, y=y_bot}
c = Vector3:new{x=x_end, y=y_bot}

d = Vector3:new{x=x_start, y=y_top}
e = Vector3:new{x=x_mid, y=y_top}
f = Vector3:new{x=x_end, y=y_top}

bot_left = Line:new{p0=a, p1=b}
bot_right = Line:new{p0=b, p1=c}
top_left = Line:new{p0=d, p1=e}
top_right = Line:new{p0=e, p1=f}

left_wall = Line:new{p0=a, p1=d}
mid_wall = Line:new{p0=b, p1=e}
right_wall = Line:new{p0=c, p1=f}

N = 400
left_grid = StructuredGrid:new{psurface=makePatch{north=top_left, east=mid_wall, south=bot_left, west=left_wall}, niv=N/2+1, njv=4}
right_grid = StructuredGrid:new{psurface=makePatch{north=top_right, east=right_wall, south=bot_right, west=mid_wall}, niv=N/2+1, njv=4}

--Create the blocks
left_block = FluidBlock:new{grid=left_grid, initialState=left, bcList={west=OutFlowBC_SimpleExtrapolate:new{xOrder=1}}}
right_block = FluidBlock:new{grid=right_grid, initialState=right, bcList={east=OutFlowBC_SimpleExtrapolate:new{xOrder=1}}}

connectBlocks(left_block, north, left_block, south)
connectBlocks(right_block, north, right_block, south)

identifyBlockConnections()

--Set timestepping
config.adjust_invalid_cell_data = true
config.max_invalid_cells = 6
config.cfl_value = 0.8
config.max_time = 0.2e-2
config.max_step = 1e7
config.dt_init = 1e-8
--config.dt_max = 0.1e-6
config.dt_plot = config.max_time / 20
