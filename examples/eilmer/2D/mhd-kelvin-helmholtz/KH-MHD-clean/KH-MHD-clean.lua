--Re-create Kelvin-Helmholtz MHD case in Eilmer4
--Initial Eilmer3 case built by Daryl Bond
--Lachlan Whyborn 9/02/2016

--set attributes of the simulation
config.title = "Kelvin-Helmholtz-MHD-eilmer4"
print(config.title)
config.dimensions = 2
config.spatial_deriv_calc = least_squares
config.flux_calculator = "hlle"
config.MHD = true
config.divergence_cleaning = true
config.axisymmetric = false
config.stringent_cfl = true

--set gas model to ideal air
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')

--create user-defined flow function
function UDF_flow(x,y,z)

	--set properties of the flow    
	rho0 = 1.2
	p0 = 101.325e3
	B0 = math.sqrt(p0)
	R = 287.0
	gam = 7.0/5.0

	L_inf = 1.0
	u_inf = math.sqrt(p0/rho0)
	t_inf = L_inf/u_inf

	rho1 = 1.0*rho0
	p1 = 50.0*p0
	T1 = p1/(rho1*R)
    
    velx1 = 5.0*(math.tanh(20.*(y+0.5))-(math.tanh(20.*(y-0.5))+1))*u_inf
    vely1 = 0.25*math.sin(2.*math.pi*x)*(math.exp(-100.*(y+0.5)*(y+0.5))-math.exp(-100.*(y-0.5)*(y-0.5)))*u_inf
    velz1 = 0.0
    
    Bx1 = 1.0*B0
    By1 = 0.0
    Bz1 = 0.0
    
    return FlowState:new{p=p1, T=T1, velx = velx1, vely = vely1, Bx = Bx1}
end
    
--==============================================================================
-- RUN DEFINITION
--==============================================================================

	rho0 = 1.2
	p0 = 101.325e3
	B0 = math.sqrt(p0)
	R = 287.0
	gam = 7.0/5.0

	L_inf = 1.0
	u_inf = math.sqrt(p0/rho0)
	t_inf = L_inf/u_inf

config.max_time = t_inf
config.max_step = 10000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/100.0

--==============================================================================
-- DOMAIN DEFINITION
--==============================================================================

start_x = 0.0
mid_x = 0.5
stop_x = 1.0

start_y = -1.0
mid_y = 0.0
stop_y = 1.0

-- nodes
a = Vector3:new{x = start_x, y = start_y}
b = Vector3:new{x = mid_x, y = start_y}
c = Vector3:new{x = stop_x, y = start_y}

d = Vector3:new{x = start_x, y = mid_y}
e = Vector3:new{x = mid_x, y = mid_y}
f = Vector3:new{x = stop_x, y = mid_y}

g = Vector3:new{x = start_x, y = stop_y}
h = Vector3:new{x = mid_x, y = stop_y}
i = Vector3:new{x = stop_x, y = stop_y}

-- lines
ab = Line:new{p0=a, p1=b}
bc = Line:new{p0=b, p1=c}
de = Line:new{p0=d, p1=e}
ef = Line:new{p0=e, p1=f}
gh = Line:new{p0=g, p1=h}
hi = Line:new{p0=h, p1=i}

ad = Line:new{p0=a, p1=d}
dg = Line:new{p0=d, p1=g}
be = Line:new{p0=b, p1=e}
eh = Line:new{p0=e, p1=h}
cf = Line:new{p0=c, p1=f}
fi = Line:new{p0=f, p1=i}

config.divB_damping_length = 1.0

--Define the blocks, boundary conditions and set the discretisation.
N = 64
grid0 = StructuredGrid:new{psurface = makePatch{north=de, east=be, south=ab, west=ad}, niv=N/2 + 1, njv=N + 1}
grid1 = StructuredGrid:new{psurface = makePatch{north=ef, east=cf, south=bc, west=be}, niv=N/2 + 1, njv=N + 1}
grid2 = StructuredGrid:new{psurface = makePatch{north=hi, east=fi, south=ef, west=eh}, niv=N/2 + 1, njv=N + 1}
grid3 = StructuredGrid:new{psurface = makePatch{north=gh, east=eh, south=de, west=dg}, niv=N/2 + 1, njv=N + 1}

blk_0 = FluidBlock:new{grid = grid0, initialState=UDF_flow, label="flow"}

blk_1 = FluidBlock:new{grid = grid1, initialState=UDF_flow, label="flow"}
                
blk_2 = FluidBlock:new{grid = grid2, initialState=UDF_flow, label="flow"}

blk_3 = FluidBlock:new{grid = grid3, initialState=UDF_flow, label="flow"}

identifyBlockConnections()

--Set periodic boundary conditions
connectBlocks(blk_0, south, blk_3, north, 0)
connectBlocks(blk_1, south, blk_2, north, 0)
connectBlocks(blk_0, west, blk_1, east, 0)
connectBlocks(blk_3, west, blk_2, east, 0)
