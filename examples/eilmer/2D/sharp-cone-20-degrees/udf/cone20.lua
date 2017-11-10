-- cone20.lua
-- Simple job-specification file for e4prep -- for use with Eilmer4
-- PJ & RG
-- 2015-02-24 -- adapted from the Python version of cone20

job_title = "Mach 1.5 flow over a 20 degree cone."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title
config.axisymmetric = true

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0, vely=0.0}

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.2, y=0.0}
c = Vector3:new{x=1.0, y=0.29118}
d = Vector3:new{x=1.0, y=1.0}
e = Vector3:new{x=0.2, y=1.0}
f = Vector3:new{x=0.0, y=1.0}
-- lower boundary including cone surface
ab = Line:new{p0=a, p1=b}; bc = Line:new{p0=b, p1=c}
-- upper boundary
fe = Line:new{p0=f, p1=e}; ed = Line:new{p0=e, p1=d}
-- vertical lines
af = Line:new{p0=a, p1=f}; be = Line:new{p0=b, p1=e}; cd = Line:new{p0=c, p1=d}
quad0 = makePatch{north=fe, east=be, south=ab, west=af}
quad1 = makePatch{north=ed, east=cd, south=bc, west=be, gridType="ao"}
-- Mesh the patches, with particular discretisation.
nx0 = 10; nx1 = 30; ny = 40
grid0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1}
-- Define the flow-solution blocks.
blk0 = FluidBlock:new{grid=grid0, initialState=inflow, label="BLOCK-0"}
blk1 = FluidBlock:new{grid=grid1, initialState=initial, label="BLOCK-1",
		      xforceList={0,0,1,0}}
setHistoryPoint{ib=1, i=9, j=0}
-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = UserDefinedBC:new{fileName="udf-bc.lua"}
blk1.bcList[east] = UserDefinedBC:new{fileName="udf-bc.lua"}
config.apply_bcs_in_parallel = false

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
-- config.dt_max = 10.0e-6
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5
