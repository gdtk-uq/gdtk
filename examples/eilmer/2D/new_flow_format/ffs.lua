-- ffs.lua
-- Forward-facing-step
-- PJ & RG  2015-03-08, 2019-07-15 simplify
-- DB 2021-04-20 test of new flow format

config.title = "Forward-facing step with supersonic flow."
print(config.title)
config.dimensions = 2


config.new_flow_format = true
config.flow_format = "eilmer4text"
-- config.flow_format = "eilmer4binary"

-- config.new_flow_format = false
-- config.flow_format = "gziptext"

config.do_temporal_DFT = true
config.do_flow_average = true


-- Gas model and flow conditions.
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=101.325e3, T=300.0}
sspeed = initial.a
print("sound speed=", sspeed)
inflow = FlowState:new{p=101.325e3, T=300.0, velx=3.0*sspeed, vely=0.0}

-- Geometry of the flow domain.
a0 = {x=0.0, y=0.0}; a1 = {x=0.0, y=0.2}; a2 = {x=0.0, y=1.0}
b0 = {x=0.6, y=0.0}; b1 = {x=0.6, y=0.2}; b2 = {x=0.6, y=1.0}
c1 = {x=3.0, y=0.2}; c2 = {x=3.0, y=1.0}
surf0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
surf1 = CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2}
surf2 = CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2}

-- Mesh the patches, with particular discretisation.
dx = 5.0e-2
nab = math.floor(0.6/dx); nbc = math.floor(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = math.floor(0.2/dx); n12 = math.floor(0.8/dx)
print("n01=", n01, "n12=", n12)
grid0 = StructuredGrid:new{psurface=surf0, niv=nab+1, njv=n01+1}
grid1 = StructuredGrid:new{psurface=surf1, niv=nab+1, njv=n12+1}
grid2 = StructuredGrid:new{psurface=surf2, niv=nbc+1, njv=n12+1}

-- Set boundary conditions that we care about.
bcList01 = { west=InFlowBC_Supersonic:new{flowState=inflow, label="inflow"} }
bcList2 = { east=OutFlowBC_Simple:new{label="outflow"} }

-- Define the flow-solution blocks and stitch them together.
blk0 = FBArray:new{grid=grid0, nib=1, njb=1, initialState=inflow, bcList=bcList01}
blk1 = FBArray:new{grid=grid1, nib=1, njb=4, initialState=inflow, bcList=bcList01}
blk2 = FBArray:new{grid=grid2, nib=4, njb=4, initialState=inflow, bcList=bcList2}
identifyBlockConnections()

-- The number of MPI tasks needs to match the mpirun command.
mpiDistributeBlocks{ntasks=3}
config.max_time = 5.0e-3  -- seconds
config.max_step = 6000
config.dt_plot = 1.0e-3
