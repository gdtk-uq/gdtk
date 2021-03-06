-- plate.lua
-- An example of using the heat transfer solver in stand-alone mode.
-- This simulation is a simple rectangular plate.
-- RG & PJ
-- 2015-04-28

job_title = "Heat transfer in a solid plate"
print(job_title)

config.dimensions = 2
config.title = job_title

-- Set dummy gas and flow properties
-- We need to fill a flow block even if we aren't interested
-- in the values.
setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}

-- Set geometry: a rectangle in the (x,y)-plane
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.1, y=0.0}
c = Vector3:new{x=0.0, y=0.02}
d = Vector3:new{x=0.1, y=0.02}
ab = Line:new{p0=a, p1=b}
ac = Line:new{p0=a, p1=c}
bd = Line:new{p0=b, p1=d}
cd = Line:new{p0=c, p1=d}
patch0 = makePatch{north=cd, east=bd, south=ab, west=ac}

-- Create a grid for the plate
nx = 100; ny = 20
grid0 = StructuredGrid:new{psurface=patch0, niv=nx+1, njv=ny+1}

-- Set a dummy flow block (note that fluid blocks need to be defined BEFORE solid blocks).
dummy = FluidBlock:new{grid=grid0, initialState=initial, label="dummy", active=false}

-- Define the solution block
-- Physical properties are those of copper
blk0 = SolidBlock:new{grid=grid0, initTemperature=300.0,
		      properties={rho=8960, k=401, Cp=386}}

-- Set boundary conditions
blk0.bcList['north'] = SolidUserDefinedBC:new{fileName='fixed-T-bc.lua'}
blk0.bcList['east'] = SolidUserDefinedBC:new{fileName='fixed-T-bc.lua'}
blk0.bcList['south'] = SolidUserDefinedBC:new{fileName='fixed-T-bc.lua'}
blk0.bcList['west'] = SolidUserDefinedBC:new{fileName='fixed-T-bc.lua'}

-- Set some simulation parameters
config.solid_domain_update_scheme = "pc"
config.max_time = 1.0 -- seconds
config.max_step = 1000
config.dt_init = 1.0e-3
config.fixed_time_step = true
config.dt_plot = 0.1


