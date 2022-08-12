-- t4m4c.lua : T4 Mach 4 nozzle steady-state variant
-- Kyle Damm, 30-07-2020

config.title = "T4 Mach4 nozzle with air in chemical equilibrium."
print(config.title)
config.dimensions = 2
config.axisymmetric = true

-- Set up inflow and initial conition.
nsp, nmodes, gm = setGasModel('cea-lut-air.lua')

-- Get flow conditions at the nozzle throat by performing ESTCj-like calculations,
-- using the CEA-backed gas model.
print("First, compute the flow conditions at the nozzle throat.")
print("shock-tube fill conditions")
gm = GasModel:new{'cea-air13species-gas-model.lua'}
state1 = GasState:new{gm}
state1.p = 131.0e3; state1.T = 300.0
gm:updateThermoFromPT(state1); gm:updateTransCoeffs(state1)
print("state1:"); printValues(state1)
--
print("normal shock, given shock speed")
Vs = 1644.0
state2, V2, Vg = gasflow.normal_shock(state1, Vs)
gm:updateThermoFromPT(state2); gm:updateTransCoeffs(state2)
print("V2=", V2, "Vg=", Vg)
print("state2:"); printValues(state2)
--
print("reflected shock")
state5, Vr = gasflow.reflected_shock(state2, Vg)
gm:updateThermoFromPT(state5); gm:updateTransCoeffs(state5)
print("Vr=", Vr)
print("state5:"); printValues(state5)
--
print("Expand from stagnation (with ratio of pressure to match observation)")
state5s, V5s = gasflow.expand_from_stagnation(state5, 8.32/23.61)
gm:updateThermoFromPT(state5s); gm:updateTransCoeffs(state5s)
print("V5s=", V5s, " Mach=", V5s/state5s.a)
print("state5s:"); printValues(state5s)
print("(h5s-h1)=", gm:enthalpy(state5s) - gm:enthalpy(state1))
--
print("Expand to throat condition (Mach 1.0001)")
state6, V6 = gasflow.expand_to_mach(state5s, 1.0001)
gm:updateThermoFromPT(state6); gm:updateTransCoeffs(state6)
print("V6=", V6, " Mach=", V6/state6.a)
print("state6:"); printValues(state6)
--
print("Quasi-one-dimensional expansion to estimated test-flow condition.")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 27.0)
gm:updateThermoFromPT(state7); gm:updateTransCoeffs(state7)
print("V7=", V7, " Mach=", V7/state7.a)
print("state7:"); printValues(state7)
--
print("Second, set up the simulation grid and the initial flow state.")

-- Set up flow state.
-- For the simulation, we are going to use the tabulated gas behaviour
-- because it will be much faster than the CEA-backed gas model.
nsp, nmodes, gm = setGasModel('cea-lut-air.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)
-- The inflow for the simulation is at the nozzle throat.
inflow = FlowState:new{p=state6.p, T=state6.T, velx=V6}
print("Species mass fractions at throat, according to CEA gas model")
for k, v in pairs(state6.ceaSavedData.massf) do
   print(string.format("massf[%s] = %g", k, v))
end

-- set the initial state to be the estimated nozzle outflow
initial = FlowState:new{p=state7.p, T=state7.T, velx=V7}

-- Define the flow domain using an imported grid.
nblocks=4
grids = {}
for i=0,nblocks-1 do
   fileName = string.format("su2-grid/block_%d_t4_nozzle.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1}
end

T_wall = 300.0 -- K
my_bcDict = {SYMMETRY = WallBC_WithSlip:new{},
             INFLOW = InFlowBC_Supersonic:new{flowCondition=inflow},
             OUTFLOW = OutFlowBC_Simple:new{},
             WALL = WallBC_NoSlip_FixedT:new{Twall=T_wall},
             METIS_INTERIOR   = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, list_mapped_cells=false}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], fillCondition=initial, bcDict=my_bcDict}
end

-- general settings
config.print_count = 20

-- invsicid flux settings
config.flux_calculator = "ausmdv"
config.extrema_clipping = false
config.interpolation_order = 2
config.thermo_interpolator = "rhop"
config.unstructured_limiter = "venkat"
config.smooth_limiter_coeff = 0.1
config.freeze_limiter_on_step = 2300

-- viscous flux settings
config.viscous = true
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.diffuse_wall_bcs_on_init = true
config.number_init_passes = 50

config.with_local_time_stepping = true
config.stringent_cfl = true

SteadyStateSolver{ 
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 100;
   start_preconditioning = 1,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   
   number_pre_steps = 10,
   number_total_steps = 1e6,
   stop_on_relative_global_residual = 1.0e-12,
   
   -- Settings for FGMRES iterative solver
   max_outer_iterations = 10,
   max_restarts = 0,

   residual_based_cfl_scheduling = true,

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 0.5,

   -- Settings for inexact Newton phase
   sigma1 = 1.0e-30,
   cfl1 = 0.5,
   tau1 = 1.0,
   eta1 = 1e-2,
   eta_strategy = "constant",
   p1 = 0.75,
   
   -- Settings control write-out
   snapshots_count = 250,
   number_total_snapshots = 5,
   write_diagnostics_count = 20,
   write_loads_count = 1000,
}
