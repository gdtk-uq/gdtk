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
config.dimensions = 3

file = io.open("case.txt", "r")
case = tonumber(file:read("*line"))
fluxCalc = tostring(file:read("*line"))
derivCalc = tostring(file:read("*line"))
xOrder = tonumber(file:read("*line"))
blocking = tostring(file:read("*line"))
ncells = tonumber(file:read("*line"))
file:close()

setGasModel('very-viscous-air.lua')
R = 287.0
p0 = 1.0e5; T0 = p0/R;
if case == 1 or case == 3 then
   u0 = 800.0; v0 = 800.0; w0 = 800.0 --u0 = 800.0; v0 = 800.0; w0 = 800.0
else
   u0 = 70.0; v0 = 90.0; w0 = 80.0
end
--initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0, velz = w0}
dofile('fill-fn.lua')

p000 = Vector3:new{x=0.0, y=0.0, z=0.0}
p100 = Vector3:new{x=1.0, y=0.0, z=0.0}
p010 = Vector3:new{x=0.0, y=1.0, z=0.0}
p110 = Vector3:new{x=1.0, y=1.0, z=0.0}
p001 = Vector3:new{x=0.0, y=0.0, z=1.0}
p101 = Vector3:new{x=1.0, y=0.0, z=1.0}
p011 = Vector3:new{x=0.0, y=1.0, z=1.0}
p111 = Vector3:new{x=1.0, y=1.0, z=1.0}
nicell = ncells; njcell = ncells; nkcell = ncells
grid = StructuredGrid:new{pvolume=TFIVolume:new{vertices={p000, p100, p110, p010, p001, p101, p111, p011}}, niv=ncells+1, njv=ncells+1, nkv=ncells+1}

bcList = {}
if case == 1 or case == 3 then
   bcList[north] = UserDefinedBC:new{fileName='udf-bc.lua'} --OutFlowBC_Simple:new{}
   bcList[east] = UserDefinedBC:new{fileName='udf-bc.lua'} --OutFlowBC_Simple:new{}
   bcList[south] = UserDefinedBC:new{fileName='udf-bc.lua'}
   bcList[west] = UserDefinedBC:new{fileName='udf-bc.lua'}
   bcList[top] = UserDefinedBC:new{fileName='udf-bc.lua'} --OutFlowBC_Simple:new{}
   bcList[bottom] = UserDefinedBC:new{fileName='udf-bc.lua'}

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
   bcList[top] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new()
      }
   }
   bcList[bottom] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new()
      }
   }
end
config.apply_bcs_in_parallel = false
if blocking == 'single' then
   blk = FluidBlock:new{grid=grid, fillCondition=gasFillFn, bcList=bcList, label="blk"}
   SBlock2UBlock(blocks[1])
else 
   blks = FluidBlockArray{grid=grid, fillCondition=gasFillFn, bcList=bcList, 
			  nib=2, njb=2, nkb=2, label="blk"}
   for i=1,8 do
      SBlock2UBlock(blocks[i])
   end 
end

--config.adjust_invalid_cell_data = true
--config.max_invalid_cells = 10
config.interpolation_order = xOrder
config.gasdynamic_update_scheme = "predictor-corrector"
config.flux_calculator = fluxCalc
config.spatial_deriv_calc = derivCalc
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
if case == 1 or case == 3 then
   config.dt_init = 1.0e-6
   config.max_time = 60.0e-3
else
   config.viscous = true
   config.dt_init = 0.5e-7
   config.max_time = 150.0e-3
   config.viscous_signal_factor = 0.1
end
config.dt_plot = config.max_time/20.0
config.max_step = 100
config.cfl_value = 0.4
config.stringent_cfl = true
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false
