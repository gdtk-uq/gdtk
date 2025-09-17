-- fp.lua : Laminar flow over a flat plate
-- Nick. N. Gibbons, based on 
--  "Implicit Large-Eddy Simulation of a Supersonic Turbulent Boundary Layer: Code Comparison"
--  Jonathan Poggie, Nicholas J. Bisek, Timothy Leger, and Ricky Tang
--  AIAA Paper, 2014

config.solver_mode = "steady"
config.dimensions = 2
config.flux_calculator = "ldfss2"
config.viscous = true

config.extrema_clipping=false
config.scale_species_after_reconstruction = false
config.enforce_species_density_positivity = false
config.epsilon_van_albada = 1e-3
config.thermo_interpolator = "rhop"

-- Gas model and flow conditions to match Table 1
nsp, nmodes, gm = setGasModel('ideal-air.gas')
delta0 = 5.375e-3 -- This is the boundary layer size further downstream, our plate is shorter
u_inf = 604.5  -- m/s
p_inf = 2.303e3  -- Pa
T_inf = 108.1  -- K
Twall = 269.5
massf_inf = {air=1.0}

dywall = 5.0e-4*delta0
L = 55*delta0 -- Truncated domain so this test runs faster
H = 0.8*delta0
H2 = 15*delta0

Nx = 750
Ny = 48

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf; gas_inf.massf = massf_inf
gm:updateThermoFromPT(gas_inf); gm:updateSoundSpeed(gas_inf); gm:updateTransCoeffs(gas_inf)

-- Set up flow state
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, massf=massf_inf}

--          outflow
--        d---------c
-- flow=> |         | -> Outflow
--        a---------b
--            wall
a = Vector3:new{x=0.0, y=0.0}; b = Vector3:new{x=L, y=0.0}
c = Vector3:new{x=L, y=H};     d = Vector3:new{x=0.0, y=H}

patch = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d}

niv = Nx+1; njv = Ny+1;
acf = dywall/H

-- Cluster towards the south (wall) boundary
cflist = {west=GeometricFunction:new{a=acf, r=1.3, N=njv},
          east=GeometricFunction:new{a=acf, r=1.3, N=njv},
}


grd = StructuredGrid:new{psurface=patch, niv=niv, njv=njv, cfList=cflist}
grid0 = registerFluidGridArray{
   grid = grd,
   nib = 4,
   njb = 1,
   fsTag="inflow",
   bcTags={east="outflow",south="wall",west="inflow"}
}


a = Vector3:new{x=0.0, y=H}; b = Vector3:new{x=L, y=H}
c = Vector3:new{x=L, y=H2};     d = Vector3:new{x=0.0, y=H2}

patch1 = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d}

cflist = {east=RobertsFunction:new{end0=true, end1=false, beta=1.05},
          west=RobertsFunction:new{end0=true, end1=false, beta=1.05},
}

grd1= StructuredGrid:new{psurface=patch1, niv=niv, njv=njv, cfList=cflist}
grid1 = registerFluidGridArray{
   grid = grd1,
   nib = 4,
   njb = 1,
   fsTag="inflow",
   bcTags={north="inflow",east="outflow",west="inflow"}
}

identifyGridConnections()
flowDict = {}
flowDict["inflow"] = inflow

bcDict = {
   wall = WallBC_NoSlip_FixedT:new{Twall=Twall, group="wall"},
   inflow = InFlowBC_Supersonic:new{flowState=inflow},
   outflow = OutFlowBC_Simple:new{}
}

makeFluidBlocks(bcDict, flowDict)

-- loads settings
config.boundary_groups_for_loads = "wall"
config.write_loads = true

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 1200,
   stop_on_relative_residual = 1.0e-8,
   number_of_phases = 1,
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = true,
   max_linear_solver_iterations = 60,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-30,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-30,
   preconditioner = "ilu",
   ilu_fill = 0,
   write_loads = true,
   total_snapshots = 2,
   steps_between_snapshots = 100,
   steps_between_diagnostics = 1,
   steps_between_loads_update = 1000000,
   write_residual_values = true,
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   frozen_preconditioner = true,
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 10,
   linear_solve_tolerance = 0.01,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.1,
   start_cfl = 2.0,
   max_cfl = 1.0e4,
   auto_cfl_exponent = 0.8
}

