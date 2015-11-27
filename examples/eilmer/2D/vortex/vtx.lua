-- vtx.lua
-- PJ & RG
-- 2015-11-27 -- adapted from the Python version of 2D/vortex

job_title = "Inviscid supersonic vortex -- flow in a bend."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
-- The following flow condition is not really important because
-- the actual data will be taken from the user-defined boundaries.
initial = FlowState:new{p=1000.0, T=348.43, velx=0.0, vely=0.0}

-- Geometry
R_inner = 1.0
R_outer = 1.384
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, R_inner}
c = Vector3:new{0.0, R_outer}
d = Vector3:new{R_inner, 0.0}
e = Vector3:new{R_outer, 0.0}
north0 = Arc:new{c, e, a}
east0 = Line:new{d, e}
south0 = Arc:new{b, d, a}
west0 = Line:new{b, c}

-- Mesh the patches, with particular discretisation.
nx = 80; ny = 40
grid0 = StructuredGrid:new{psurface=makePatch{north0, east0, south0, west0},
			   niv=nx+1, njv=ny+1}
-- Define the flow-solution block.
bcList = {north=UserDefinedBC:new{fileName='udf-vortex-flow.lua'},
	  east=OutFlowBC_Simple:new{},
	  south=UserDefinedBC:new{fileName='udf-vortex-flow.lua'},
	  west=UserDefinedBC:new{fileName='udf-vortex-flow.lua'}}
blk0 = SBlock:new{grid=grid0, fillCondition=initial, bcList=bcList, label="Duct"}

-- Do a little more setting of global data.
config.max_time = 20.0e-3  -- seconds
config.max_step = 6000
config.dt_init = 1.0e-6
config.dt_plot = 5.0e-3
