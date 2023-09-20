-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid
--
--          2015-06-03
--          Added scripting to configure case based on 'case.txt'
--
--          2023-08-21, RJG and KAD
--          Updated for use with lmr5, specifically:
--          1. staged prep;
--          2. new steady-state solver input
--          3. case configuration via 'config.txt'
--

config.title = "Method of Manufactured Solutions."
print(config.title)
config.dimensions = 3

-- Case is configured by a higher-level controller.
-- Options are placed in config.txt as simple
-- variable assignment statements and then used
-- in this job script.
dofile("config.txt")

-- Pull gas constants from constants.txt
dofile("constants.txt")

setGasModel('very-viscous-air.lua')
rho0 = 1.0; p0 = 1.0e5; T0 = p0/(rho0*Rgas);
u0 = 70.0; v0 = 90.0; w0 = 80.0
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0, velz=w0}
flowStates = {initial=initial}

p000 = Vector3:new{x=0.0, y=0.0, z=0.0}
p100 = Vector3:new{x=1.0, y=0.0, z=0.0}
p010 = Vector3:new{x=0.0, y=1.0, z=0.0}
p110 = Vector3:new{x=1.0, y=1.0, z=0.0}
p001 = Vector3:new{x=0.0, y=0.0, z=1.0}
p101 = Vector3:new{x=1.0, y=0.0, z=1.0}
p011 = Vector3:new{x=0.0, y=1.0, z=1.0}
p111 = Vector3:new{x=1.0, y=1.0, z=1.0}
nicell = ncells; njcell = ncells; nkcell = ncells
sgrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices={p000, p100, p110, p010, p001, p101, p111, p011}},
                           niv=ncells+1, njv=ncells+1, nkv=ncells+1}

-- All faces have UDF, so just loop over all structured face names and assign UDF
bcTags = {}
for face,idx in pairs(Face) do
   bcTags[idx] = "udf"
end

grid = registerGrid{
   grid=UnstructuredGrid:new{sgrid=sgrid},
   fsTag="initial",
   bcTags=bcTags
}

udfBC = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivActionAtBndryFaces = { UserDefinedInterface:new{fileName='udf-bc.lua'},
                                         UpdateThermoTransCoeffs:new() }
}
bcs = {udf=udfBC}
config.apply_bcs_in_parallel = false

makeFluidBlocks(bcs, flowStates)

config.include_ghost_cells_in_spatial_deriv_clouds = true
config.interpolation_order = interpolation_order
config.flux_calculator = flux_calculator
config.spatial_deriv_calc = spatial_deriv_calc
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.viscous = true
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 70,
   stop_on_relative_residual = 1e-15,
   number_of_phases = 2,
   phase_changes_at_steps = { 20 },
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = false,
   max_linear_solver_iterations = 50,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 2,
   steps_between_status = 5,
   steps_between_snapshots = 20,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   frozen_preconditioner = true,
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 10,
   linear_solve_tolerance = 0.1,
   use_local_timestep = true,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.9,
   start_cfl = 100.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 100.0,
}

