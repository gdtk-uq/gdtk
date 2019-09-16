-- vtx.lua
-- PJ & RG
-- Exercises the user-defined boundary condition. 
-- 2015-11-27 adapted from the Python version of 2D/vortex
-- 2019-07-15 refresh

config.title = "Inviscid supersonic vortex -- flow in a bend."
print(config.title)
config.dimensions = 2

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=1000.0, T=348.43, velx=0.0, vely=0.0}

-- Geometry
R_inner = 1.0
R_outer = 1.384
a = {x=0.0, y=0.0}
b = {x=0.0, y=R_inner}
c = {x=0.0, y=R_outer}
d = {x=R_inner, y=0.0}
e = {x=R_outer, y=0.0}
patch = makePatch{north=Arc:new{p0=c, p1=e, centre=a},
		  east=Line:new{p0=d, p1=e},
		  south=Arc:new{p0=b, p1=d, centre=a},
		  west=Line:new{p0=b, p1=c}}
nx = 80; ny = 40
grid0 = StructuredGrid:new{psurface=patch, niv=nx+1, njv=ny+1}
-- Flow domain
--[[
bcList = {north=UserDefinedBC:new{fileName='udf-vortex-flow.lua'},
	  east=OutFlowBC_Simple:new{},
	  south=UserDefinedBC:new{fileName='udf-vortex-flow.lua'},
	  west=UserDefinedBC:new{fileName='udf-vortex-flow.lua'}}
--]]
bcList = {north=WallBC_WithSlip1:new{},
	  east=OutFlowBC_Simple:new{},
	  south=WallBC_WithSlip1:new{},
	  west=UserDefinedBC:new{fileName='udf-vortex-flow.lua'}}
blk = FluidBlockArray{grid=grid0, nib=4, njb=1,
                      initialState=initial, bcList=bcList,
		      label="Duct"}

config.max_time = 20.0e-3
config.max_step = 6000
config.dt_init = 1.0e-6
config.dt_plot = 5.0e-3
