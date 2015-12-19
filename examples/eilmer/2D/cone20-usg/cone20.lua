-- cone20.lua
-- Unstructured Grid Example -- for use with Eilmer4
-- PJ & RG
-- 2015-11-08 -- adapted from cone20-simple

job_title = "Mach 1.5 flow over a 20 degree cone -- Unstructured Grid."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title
config.axisymmetric = true

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0, vely=0.0}

-- Demo: Verify Mach number of inflow.
Q = GasState:new{gm}
Q.p = 95.84e3
Q.T = {1103.0}
print("T", Q.T[1])
Q.massf = {air=1.0}
gm:updateSoundSpeed(Q)
print("Sound speed= ", Q.a)
print("Inflow Mach number= ", 1000.0/Q.a)

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.2, 0.0}
c = Vector3:new{1.0, 0.29118}
d = Vector3:new{1.0, 1.0}
e = Vector3:new{0.2, 1.0}
f = Vector3:new{0.0, 1.0}
ab = Line:new{a, b}; bc = Line:new{b, c} -- lower boundary including cone surface
fe = Line:new{f, e}; ed = Line:new{e, d} -- upper boundary
af = Line:new{a, f}; be = Line:new{b, e}; cd = Line:new{c, d} -- vertical lines
-- Mesh the patches, with particular discretisation.
nx0 = 10; nx1 = 30; ny = 40
grid0 = StructuredGrid:new{psurface=makePatch{fe, be, ab, af}, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=makePatch{ed, cd, bc, be, gridType="ao"}, niv=nx1+1, njv=ny+1}
-- Define the flow-solution blocks.
blk0 = UBlock:new{grid=UnstructuredGrid:new{sgrid=grid0},
		  fillCondition=inflow, label="BLOCK-0"}
blk1 = UBlock:new{grid=UnstructuredGrid:new{sgrid=grid1}, 
		  fillCondition=initial,
		  label="BLOCK-1",
		  hcellList={9,0}}
-- Set boundary conditions.
-- identifyBlockConnections()
blk0.bcList[3] = InFlowBC_Supersonic:new{flowCondition=inflow, label="inflow-boundary"}
blk0.bcList[1] = ExchangeBC_MappedCell:new{}
blk1.bcList[3] = ExchangeBC_MappedCell:new{}
blk1.bcList[1] = OutFlowBC_Simple:new{label="outflow-boundary"}

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
-- config.dt_max = 10.0e-6
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5
