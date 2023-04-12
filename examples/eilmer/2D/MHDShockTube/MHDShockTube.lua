--MHD Shock tube example
-- Conditions from Dai and Woodward, the
-- third shock tube problem in
-- "A Simple Finite Difference Scheme for
-- Multidimensional Magnetohydrodynamical
-- Equations", JCP 142. Figure 7 and Table 2.

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
nsp, nmodes = setGasModel('ideal-He.lua')
gmodel = GasModel:new{'ideal-He.lua'}

-- From Dai and Woodward
Bx = 5 / math.sqrt(4 * math.pi)

rho_left = 1
p_left = 1
vx_left = 11.7106
vy_left = -0.312542
vz_left = -0.495260
By_left = 1 / math.sqrt(4 * math.pi)
Bz_left = -1 / math.sqrt(4 * math.pi)
Q_left = GasState:new{gmodel}
Q_left.p = p_left; Q_left.rho = rho_left
gmodel:updateThermoFromRHOP(Q_left)
T_left = Q_left.T
left = FlowState:new{T=T_left, p=p_left, velx=vx_left, vely=vy_left, velz=vz_left, Bx=Bx, By=By_left, Bz=Bz_left}

rho_right = 1
p_right = 1
vx_right = -11.7106
vy_right = 0.312542
vz_right = -0.495260
By_right = 1 / math.sqrt(4 * math.pi)
Bz_right = 1 / math.sqrt(4 * math.pi)
Q_right = GasState:new{gmodel}
Q_right.p = p_right; Q_right.rho = rho_right
gmodel:updateThermoFromRHOP(Q_right)
T_right = Q_right.T
right = FlowState:new{T=T_right, p=p_right, velx=vx_right, vely=vy_right, velz=vz_right, Bx=Bx, By=By_right, Bz=Bz_right}

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

N = 512
left_grid = StructuredGrid:new{psurface=makePatch{north=top_left, east=mid_wall, south=bot_left, west=left_wall}, niv=N/2+1, njv=4}
right_grid = StructuredGrid:new{psurface=makePatch{north=top_right, east=right_wall, south=bot_right, west=mid_wall}, niv=N/2+1, njv=4}

--Create the blocks
left_block = FluidBlock:new{grid=left_grid, initialState=left, bcList={west=OutFlowBC_SimpleExtrapolate:new{xOrder=1}}}
right_block = FluidBlock:new{grid=right_grid, initialState=right, bcList={east=OutFlowBC_SimpleExtrapolate:new{xOrder=1}}}

connectBlocks(left_block, "north", left_block, "south")
connectBlocks(right_block, "north", right_block, "south")

identifyBlockConnections()

--Set timestepping
config.adjust_invalid_cell_data = false
config.max_invalid_cells = 0
config.cfl_value = 0.5
config.max_time = 6e-2
config.max_step = 1e6
config.dt_init = 1e-5
--config.dt_max = 0.1e-6
config.dt_plot = config.max_time / 5
