-- odw.lua
-- Powers and Aslam oblique detonation wave verification case.
-- PJ & RG
-- 2017-01-07: ported from the Eilmer3 example

-- We can set individual attributes of the global data object.
config.title = "Oblique detonation wave with Powers-Aslam gas model."
print(config.title)
config.dimensions = 2

nsp, nmodes, gm = setGasModel('powers-aslam-gas-model.lua')
print("GasModel has nsp= ", nsp, " nmodes= ", nmodes)
massf1 = {A=1.0, B=0.0}
initial = FlowState:new{p=28.7e3, T=300.0, massf=massf1}
inflow = FlowState:new{p=86.1e3, T=300.0, velx=964.302, massf=massf1}

-- Geometry
xmin = -0.25
xmax = 1.75
ymin = 0.0
ymax = 2.0

dofile("analytic.lua")
myWallFn = create_wall_function(0.0, xmax)

-- Set up two patches in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=xmin, y=0.0}
b = Vector3:new{x=0.0, y=0.0}
c = Vector3:new{x=myWallFn(1.0).x, y=myWallFn(1.0).y}
d = Vector3:new{x=xmin, y=ymax}
e = Vector3:new{x=0.0, y=ymax}
f = Vector3:new{x=xmax, y=ymax}
south0 = Line:new{p0=a, p1=b} -- upstream of wedge
south1 = LuaFnPath:new{luaFnName="myWallFn"} -- wedge surface
north0 = Line:new{p0=d, p1=e}; north1 = Line:new{p0=e, p1=f}
west0 = Line:new{p0=a, p1=d} -- inflow boundary
east0west1 = Line:new{p0=b, p1=e} -- vertical line, between patches
east1 = Line:new{p0=c, p1=f} -- outflow boundary
patch0 = makePatch{north=north0, east=east0west1, south=south0, west=west0}
patch1 = makePatch{north=north1, east=east1, south=south1, west=east0west1}
-- Mesh the patches, with particular discretisation.
factor = 2 -- for adjusting the grid resolution
nxcells = math.floor(40*factor)
nycells = math.floor(40*factor)
fraction0 = (0-xmin)/(xmax-xmin) -- fraction of domain upstream of wedge
nx0 = math.floor(fraction0*nxcells); nx1 = nxcells-nx0; ny = nycells
grid0 = StructuredGrid:new{psurface=patch0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=patch1, niv=nx1+1, njv=ny+1}
-- Define the flow-solution blocks and set boundary conditions.
-- We split the patches into roughly equal blocks so that
-- we make good use of our multicore machines.
blk0 = FBArray:new{grid=grid0, initialState=inflow, nib=1, njb=2,
                   bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                           north=OutFlowBC_Simple:new{}}}
blk1 = FBArray:new{grid=grid1, initialState=initial, nib=7, njb=2,
                   bcList={east=OutFlowBC_Simple:new{},
                           north=OutFlowBC_Simple:new{}}}
identifyBlockConnections()

-- Do a little more setting of global data.
config.reacting = true
config.reactions_file = "powers-aslam-gas-model.lua"
config.flux_calculator = "ausmdv"
config.max_time = 2.0e-2  -- seconds
config.max_step = 300000
config.dt_init = 1.0e-6
config.dt_plot = config.max_time / 40
config.dt_history = 10.0e-5

dofile("sketch-domain.lua")
