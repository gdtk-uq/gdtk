-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid
--          2015-06-03
--          Added scripting to configure case based on 'case.txt'
--          Ported to steady-state solver (Kyle A. Damm)

config.title = "Method of Manufactured Solutions for k-omega turbulence model."
print(config.title)
config.dimensions = 3

config.turbulence_model = "k_omega"

file = io.open("case.txt", "r")
fluxCalc = tostring(file:read("*line"))
derivCalc = tostring(file:read("*line"))
derivLocation = tostring(file:read("*line"))
xOrder = tonumber(file:read("*line"))
blocking = tostring(file:read("*line"))
ncells = tonumber(file:read("*line"))
file:close()

setGasModel('very-viscous-air.lua')
R = 287.0
p0 = 1.0e5; T0 = p0/R;
u0 = 70.0; v0 = 90.0; w0 = 80.0
tke0 = 780.0; omega0 = 150.0
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0, velz=w0, tke=tke0, omega=omega0}
--dofile("fill-fn.lua")


p000 = Vector3:new{x=0.0, y=0.0, z=0.0}
p100 = Vector3:new{x=1.0, y=0.0, z=0.0}
p010 = Vector3:new{x=0.0, y=1.0, z=0.0}
p110 = Vector3:new{x=1.0, y=1.0, z=0.0}
p001 = Vector3:new{x=0.0, y=0.0, z=1.0}
p101 = Vector3:new{x=1.0, y=0.0, z=1.0}
p011 = Vector3:new{x=0.0, y=1.0, z=1.0}
p111 = Vector3:new{x=1.0, y=1.0, z=1.0}
nicell = ncells; njcell = ncells; nkcell = ncells
grid = StructuredGrid:new{pvolume=TFIVolume:new{vertices={p000, p100, p110, p010, p001, p101, p111, p011}}, niv=nicell+1, njv=njcell+1, nkv = nkcell+1}

bcList = {}
bcList['north'] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
bcList['east'] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
bcList['south'] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
bcList['west'] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
bcList['top'] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
bcList['bottom'] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
config.apply_bcs_in_parallel = false
if blocking == 'single' then
   blk = FluidBlock:new{grid=grid, fillCondition=initial, bcList=bcList,
		    label='blk'}
   SBlock2UBlock(fluidBlocks[1])
else
   blks = FBArray:new{grid=grid, fillCondition=initial, bcList=bcList,
		      nib=2, njb=2, nkb=2, label="blk"}
   for i=1,8 do
      SBlock2UBlock(fluidBlocks[i])
   end
end

config.include_ghost_cells_in_spatial_deriv_clouds = true
config.interpolation_order = xOrder
config.flux_calculator = fluxCalc
config.spatial_deriv_calc = derivCalc
config.spatial_deriv_locn = "cells"
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.viscous = true
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false
-- Set simulation parameters
config.dt_loads = 1.0e-3
config.print_count = 100
SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 10,
   start_preconditioning = 1,

   use_scaling = true,
   use_complex_matvec_eval = true,

   number_pre_steps = 1,
   number_total_steps = 100000,
   stop_on_relative_global_residual = 1.0e-10,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 20,
   max_restarts = 0,

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 1.0,
   eta0 = 0.1,
   tau0 = 1.0,
   sigma0 = 1.0e-30,
   p0 = 1.0,

   -- Settings for inexact Newton phase
   cfl1 = 1.0,
   tau1 = 1.0,
   sigma1 = 1.0e-30,
   p1 = 1.0,
   eta1 = 0.1,
   eta1_min = 0.1,
   eta_strategy = "geometric",

   -- Settings control write-out
   snapshots_count = 50,
   number_total_snapshots = 5,
   write_diagnostics_count = 20
}
