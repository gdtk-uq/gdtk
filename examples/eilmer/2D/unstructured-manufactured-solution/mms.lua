-- mms.lua
--
-- Authors: Kyle Damm
-- Date: 2016-05-03
-- History: Ported from eilmer4 structured example
--          Presently only Euler case on regular unstructured grid
--
--             UDF
--        p01--------p11
--        |           |
--   UDF  |           | UDF
--        |           |
--        |           |
--        p00--------p10
--            UDF
--
--
config.title = "Method of Manufactured Solutions."
print(config.title)
config.dimensions = 2

-- test parameters

fluxCalc = "ausmdv"
derivCalc = "divergence"
xOrder = 2
blocking = "single"
ncells = 16

-- gas properties and intitial flowstate

setGasModel('very-viscous-air.lua')
R = 287.0
p0 = 1.0e5; T0 = p0/R;
u0 = 800.0; v0 = 800.0
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}

-- geometry

p00 = Vector3:new{x=0.0, y=0.0}
p10 = Vector3:new{x=1.0, y=0.0}
p01 = Vector3:new{x=0.0, y=1.0}
p11 = Vector3:new{x=1.0, y=1.0}
nicell = ncells; njcell = ncells
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
			  niv=nicell+1, njv=njcell+1}

-- BC's

bcList = {}
bcList[0] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList[1] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList[2] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList[3] = UserDefinedBC:new{fileName='udf-bc.lua'}
config.apply_bcs_in_parallel = false

-- setup blocks

if blocking == 'single' then
    blk = UBlock:new{grid=UnstructuredGrid:new{sgrid=grid}, fillCondition=initial, bcList=bcList,
		     label='blk'}
else 
   blks = UBlockArray{grid=UnstructuredGrid:new{sgrid=grid}, fillCondition=initial, bcList=bcList, 
		      nib=2, njb=2, label="blk"}
end

-- configuration parameters

config.interpolation_order = xOrder
config.gasdynamic_update_scheme = "predictor-corrector"
config.flux_calculator = fluxCalc
config.spatial_deriv_calc = derivCalc
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.dt_init = 1.0e-6
config.max_time = 60.0e-3
config.dt_plot = config.max_time/20.0
config.max_step = 3000000
config.cfl_value = 0.5
config.stringent_cfl = true
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false
