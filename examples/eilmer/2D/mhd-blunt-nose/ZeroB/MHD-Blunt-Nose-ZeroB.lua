--Simulation of the Blunt Nose MHD Experiment run by D. Gildfind
--Lachlan Whyborn 06/04/2016

--set attributes of the simulation
config.title = "MHD-Blunt-Nose"
config.dimensions = 2
config.flux_calculator = "hlle"
config.MHD = true
config.axisymmetric = false
config.divergence_cleaning = true

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')

--create initial static flow state
P0 = 1.
T0 = 220.

function Static(x, y, z)
	 --This will describe the magnetic field created by a 2D magnetic dipole, modelled by superimposing the magnetic field of 2 infinite length current carrying wires. Returns a flowstate object.
	 --Set static properties
	 P0 = 5.
	 T0 = 220.

	 --Location of wires
	 wire1x = 0.
	 wire1y = 1e-6
	 wire2x = 0.
	 wire2y = -1e-6

	 --Current in wires (scale until resultant magnetic field is at desired magnitude
	 I = 0.

	 --Magnetic Field
	 r1 = ((x - wire1x)^2 + (y - wire1y)^2)^0.5
	 r2 = ((x - wire2x)^2 + (y - wire2y)^2)^0.5
	 B01 = I / (2. * math.pi * r1)
	 B02 = -I / (2. * math.pi * r2)

	 B0y = (x - wire1x) / r1 * B01 + (x - wire2x) / r2 * B02
	 B0x = -(y - wire1y) * B01 / r1 - (y - wire2y) / r2 * B02
	 return FlowState:new{p = P0, T = T0, Bx = B0x, By = B0y}
end

--Create the inflow state
velin = 6283.
Tin = 438.
Pin = 53.

Inflow = FlowState:new{p = Pin, T = Tin, velx = velin}

--Create simulation domain- currently 2D case
--Radius of blunt nose
r = 0.0381/2
centre = Vector3:new{x=0.0, y=0.0}

--Create Nodes
a = Vector3:new{x=-r, y=0.0}
b = Vector3:new{x=0.0, y=r}
c = Vector3:new{x=r, y=0.0}
d = Vector3:new{x=-6*r, y= 0.0}
e = Vector3:new{x=0.0, y=6*r}
f = Vector3:new{x=6*r, y=0.0}

--Create lines for object surface
ab = Arc:new{p0=a, p1=b, centre=centre}
cb = Arc:new{p0=c, p1=b, centre=centre}

--Create lines for top of domain
de = Arc:new{p0=d, p1=e, centre=centre}
fe = Arc:new{p0=f, p1=e, centre=centre}

--Connecting Lines
da = Line:new{p0=d, p1=a}
eb = Line:new{p0=e, p1=b}
be = Line:new{p0=b, p1=e}
cf = Line:new{p0=c, p1=f}

--Set the damping length for divergence cleaning
config.divB_damping_length = 1. / (6 * r)

--Create clustering functions
cluster = RobertsFunction:new{end0=false, end1=true, beta=1.01}
inversecluster = RobertsFunction:new{end0=true, end1=false, beta=1.01}

--Create patches
nx = 100
ny = 100
grid1 = StructuredGrid:new{psurface = makePatch{north=eb, east=ab, south=da, west=de}, cfList = {north=cluster, south=cluster}, niv=nx + 1, njv=ny + 1}
grid2 = StructuredGrid:new{psurface = makePatch{north=be, east=fe, south=cf, west=cb}, cfList = {north=inversecluster, south=inversecluster},niv= nx + 1, njv=ny + 1}

--Create blocks
blk_1 = FluidBlock:new{grid = grid1, fillCondition = Static}
blk_2 = FluidBlock:new{grid = grid2, fillCondition = Static}

--Set boundary conditions
identifyBlockConnections()
blk_1.bcList[east] = WallBC_WithSlip:new{}
blk_1.bcList[west] = InFlowBC_Supersonic:new{flowCondition = Inflow}
blk_2.bcList[west] = WallBC_WithSlip:new{}
blk_2.bcList[east] = OutFlowBC_Simple:new{}

--set global time-stepping data
config.max_step = 1e7
config.dt_init = 1.0e-8
config.dt_max = 1.0e-5
config.max_time = 2.0e-4
config.dt_plot = 0.05 * config.max_time
