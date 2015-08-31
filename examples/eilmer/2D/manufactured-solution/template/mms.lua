-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid
--          2015-06-03
--          Added scripting to configure case based on 'case.txt'

config.title = "Method of Manufactured Solutions."
print(config.title)
config.dimensions = 2

file = io.open("case.txt", "r")
case = tonumber(file:read("*line"))
fluxCalc = tostring(file:read("*line"))
xOrder = tonumber(file:read("*line"))
blocking = tostring(file:read("*line"))
ncells = tonumber(file:read("*line"))
file:close()

setGasModel('very-viscous-air.lua')
R = 287.0
p0 = 1.0e5; T0 = p0/R;
if case == 1 or case == 3 then
   u0 = 800.0; v0 = 800.0
else
   u0 = 70.0; v0 = 90.0
end
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}

p00 = Vector3:new{0.0, 0.0}
p10 = Vector3:new{1.0, 0.0}
p01 = Vector3:new{0.0, 1.0}
p11 = Vector3:new{1.0, 1.0}
nicell = ncells; njcell = ncells
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
			  niv=nicell+1, njv=njcell+1}

bcList = {}
if case == 1 or case == 3 then
   bcList[north] = ExtrapolateOutBC:new{xOrder=1}
   bcList[east] = ExtrapolateOutBC:new{xOrder=1}
   bcList[south] = UserDefinedBC:new{fileName='udf-bc.lua'}
   bcList[west] = UserDefinedBC:new{fileName='udf-bc.lua'}
else
   bcList[north] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new()
      }
   }
   bcList[east] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new()
      }
   }
   bcList[south] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new()
      }
   }
   bcList[west] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new()
      }
   }
end
config.apply_bcs_in_parallel = false
if blocking == 'single' then
    blk = SBlock:new{grid=grid, fillCondition=initial, bcList=bcList,
		     label='blk'}
else 
   blks = SBlockArray{grid=grid, fillCondition=initial, bcList=bcList, 
		      nib=2, njb=2, label="blk"}
end

config.interpolation_order = xOrder
config.gasdynamic_update_scheme = "predictor-corrector"
config.flux_calculator = fluxCalc
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
if case == 1 or case == 3 then
   config.dt_init = 1.0e-6
   config.max_time = 60.0e-3
else
   config.viscous = true
   config.dt_init = 1.0e-7
   config.max_time = 150.0e-3
   config.viscous_signal_factor = 0.1
end
config.dt_plot = config.max_time/20.0
config.max_step = 3000000
config.cfl = 0.5
config.stringent_cfl = 1
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false

				 
		
