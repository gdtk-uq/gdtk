-- ffs.lua -- Forward-facing-step example
-- PJ & RG  2015-03-08

-- We can set individual attributes of the global data object.
config.title = "Forward-facing step with supersonic flow."
print(config.title)
config.dimensions = 2

-- Gas model and flow conditions.
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=101.325e3, T=300.0}
sspeed = initial:toTable().a
print("sound speed=", sspeed)
inflow = FlowState:new{p=101.325e3, T=300.0, velx=3.0*sspeed, vely=0.0}

-- Geometry of the flow domain.
a0 = Vector3:new{0.0,0.0}; a1 = Vector3:new{0.0,0.2}; a2 = Vector3:new{0.0,1.0}
b0 = Vector3:new{0.6,0.0}; b1 = Vector3:new{0.6,0.2}; b2 = Vector3:new{0.6,1.0}
c1 = Vector3:new{3.0,0.2}; c2 = Vector3:new{3.0,1.0}
surf0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
surf1 = CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2}
surf2 = CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2}

-- Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = math.floor(0.6/dx); nbc = math.floor(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = math.floor(0.2/dx); n12 = math.floor(0.8/dx)
print("n01=", n01, "n12=", n12)
grid0 = StructuredGrid:new{psurface=surf0, niv=nab+1, njv=n01+1}
grid1 = StructuredGrid:new{psurface=surf1, niv=nab+1, njv=n12+1}
grid2 = StructuredGrid:new{psurface=surf2, niv=nbc+1, njv=n12+1}

-- Set boundary conditions that we care about.
bcList0 = {north=nil, east=nil, south=nil, 
	   west=SupInBC:new{flowCondition=inflow, label="inflow-boundary"}}
bcList1 = {north=nil, east=nil, south=nil,
	   west=SupInBC:new{flowCondition=inflow, label="inflow-boundary"}}
bcList2 = {north=nil, east=ExtrapolateOutBC:new{label="outflow-boundary"},
	   south=nil, west=nil}

-- Define the flow-solution blocks and stitch them together.
blk0 = SBlockArray{grid=grid0, nib=1, njb=1, 
		   fillCondition=inflow, bcList=bcList0, label="BLOCK-0"}
blk1 = SBlockArray{grid=grid1, nib=1, njb=4,
		   fillCondition=inflow, bcList=bcList1, label="BLOCK-1"}
blk2 = SBlockArray{grid=grid2, nib=4, njb=4,
		   fillCondition=inflow, bcList=bcList2, label="BLOCK-2"}
identifyBlockConnections()

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 6000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
-- config.dt_max = 10.0e-6
config.dt_plot = 1.0e-3
config.dt_history = 10.0e-6
