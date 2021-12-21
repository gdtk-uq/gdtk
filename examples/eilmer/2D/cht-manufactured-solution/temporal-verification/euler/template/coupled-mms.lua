title = "Method of Manufactured Solutions, coupled domains case."
print(title)
config.dimensions = 2
config.title = title

-- Parse dt value from file
file = io.open("case.txt", "r")
dt = tonumber(file:read("*line"))
file:close()

-- read constants from file
dofile('constants.txt')

setGasModel('very-viscous-air.lua')

-- Load initial condition fill functions from
-- auto-generated files
dofile('fill-fn.lua')
dofile('fill-solid-fn.lua')


a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=L, y=0.0}
c = Vector3:new{x=0.0, y=H}
d = Vector3:new{x=L, y=H}
e = Vector3:new{x=0.0, y=1.5*H}
f = Vector3:new{x=L, y=1.5*H}
patch0 = makePatch{north=Line:new{p0=c, p1=d},
		   east=Line:new{p0=b, p1=d},
		   south=Line:new{p0=a, p1=b},
		   west=Line:new{p0=a, p1=c}}
patch1 = makePatch{north=Line:new{p0=e, p1=f},
		   east=Line:new{p0=d, p1=f},
		   south=Line:new{p0=c, p1=d},
		   west=Line:new{p0=c, p1=e}}
ncells = 16
nx0 = ncells
ny0 = ncells
ny1 = math.floor(ny0/2)

grid0 = StructuredGrid:new{psurface=patch0, niv=nx0+1, njv=ny0+1}
grid1 = StructuredGrid:new{psurface=patch1, niv=nx0+1, njv=ny1+1}

blk0 = FluidBlock:new{grid=grid0, initialState=gasFillFn, label="blk0"}
blk0.bcList['east'] = UserDefinedBC:new{fileName='udf-bc.lua'}
blk0.bcList['south'] = UserDefinedBC:new{fileName='udf-bc.lua'}
blk0.bcList['west'] = UserDefinedBC:new{fileName='udf-bc.lua'}

blk1 = SolidBlock:new{grid=grid1, initTemperature=solidFillFn,
		      properties={rho=rho_s, k=k_s, Cp=Cp_s}}
blk1.bcList['north'] = SolidUserDefinedBC:new{fileName='udf-solid-bc.lua'}
blk1.bcList['east'] = SolidUserDefinedBC:new{fileName='udf-solid-bc.lua'}
blk1.bcList['west'] = SolidUserDefinedBC:new{fileName='udf-solid-bc.lua'}

identifyBlockConnections()

config.interpolation_order = 2
config.gasdynamic_update_scheme = 'euler'
config.flux_calculator = 'ausm_plus_up'
config.M_inf = 0.1
config.viscous = true
config.spatial_deriv_calc = 'divergence'
config.max_time = 0.001 - 1.0e-8 -- to allow for round-off when trying to hit 0.001
config.max_step = 10000000
config.dt_init = dt
config.fixed_time_step = true
config.dt_plot = config.max_time/20

config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.udf_solid_source_terms = true
config.udf_solid_source_terms_file = 'udf-solid-source-terms.lua'
