-- nz.lua : T4 Mach 7 nozzle steady-state, lmr example 
-- References:
--   "Experimental validation of the T4 Mach 7.0 nozzle"
--   W. Y. K. Chan, M. K. Smart and P. A. Jacobs
--   School of Mechanical & Mining Engineering
--   Research Report Research Report Number 2014/14
--   https://espace.library.uq.edu.au/view/UQ:378568

--   "Flowpath design of an axisymmetric Mach 7.0 nozzle for T4".
--   Chan, W. Y. K., Smart, M. K., and Jacobs, P. A.
--   Technical Report, 2013/02.
--   School of Mechanical and Mining Engineering, UQ
--   https://espace.library.uq.edu.au/view/UQ:378572 
-- Nick Gibbons, 26/05/26

config.solver_mode = 'steady'
config.dimensions = 2
config.axisymmetric = true
config.viscous = true
config.flux_calculator = "ausmdv"
config.turbulence_model = "spalart_allmaras_edwards"
config.epsilon_van_albada= 1e-3
config.extrema_clipping=false

nsp, nmodes, gm = setGasModel('tp-air.lua')

-- Generate inflow conditions by the same procedure as in nenzf1d.
-- Use shot 11311 from Table 1 in Chan et al.
Tfill = 300.0
pfill = 200.0e3
Vs = 1679.0
pe = 19.33e6

-- Set the number of cores to parallelise the simulation here.
ncores = 16

state1 = GasState:new{gm}
state1.p = pfill; state1.T = Tfill
gm:updateThermoFromPT(state1); gm:updateTransCoeffs(state1)

state2, V2, Vg = gasflow.normal_shock(state1, Vs)
gm:updateThermoFromPT(state2); gm:updateTransCoeffs(state2)

state5, Vr = gasflow.reflected_shock(state2, Vg)
state5s, V5s = gasflow.expand_from_stagnation(state5, pe/state5.p)
gm:updateThermoFromPT(state5s); gm:updateTransCoeffs(state5s)

state6, V6 = gasflow.expand_to_mach(state5s, 1.0001)

state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 169.2)

print("state1:"); printValues(state1)
print("state2:"); printValues(state2); print("V: ", V2, "\n")
print("state5:"); printValues(state5); print("\n")
print("state5s:"); printValues(state5s); print("\n")
print("state6:"); printValues(state6); print("V: ", V6, "\n")
print("state7:"); printValues(state7); print("V: ", V7, "\n")

-- Set up inflow and initial conition.
nsp, nmodes, gm = setGasModel('tp-air.lua')

gas_inf = GasState:new{gm}
gas_inf.p = state6.p;
gas_inf.T = state6.T;
gm:updateThermoFromPT(gas_inf)
gm:updateTransCoeffs(gas_inf)
print("tpgm state 6: ", gas_inf.p, gas_inf.T, gas_inf.rho, gas_inf.massf.O2)

-- Estimate turbulence quantities
turb_lam_viscosity_ratio = 5.0 -- Fully turbulent, equation (8) from Allmaras (2012)
nu_inf = gas_inf.mu/gas_inf.rho
nuhat_inf = turb_lam_viscosity_ratio*nu_inf

inflow = FlowState:new{p=gas_inf.p, T=gas_inf.T, velx=V6, massf=gas_inf.massf, nuhat=nuhat_inf}
outflow = FlowState:new{p=state7.p, T=state7.T, velx=V7, massf=gas_inf.massf, nuhat=nuhat_inf}
print("inflow: ", inflow)
print("outflow: ", outflow)

function initial(x,y,z)
   if (x<0.075) then
       return inflow
   else
       return outflow
   end
end

npoints = 9
ctrl_pts = {
  Vector3:new{x= 0.000000e+00, y=1.050000e-02},
  Vector3:new{x= 3.000000e-02, y=1.050000e-02},
  Vector3:new{x= 1.200000e-01, y=6.750100e-02},
  Vector3:new{x= 1.600000e-01, y=5.080100e-02},
  Vector3:new{x= 2.600000e-01, y=9.194200e-02},
  Vector3:new{x= 4.000000e-01, y=1.103570e-01},
  Vector3:new{x= 5.700000e-01, y=1.179570e-01},
  Vector3:new{x= 7.700000e-01, y=1.355220e-01},
  Vector3:new{x= 1.000000e+00, y=1.365830e-01},
}

ctrl_pts2 = {
  Vector3:new{x= 0.000000e+00, y=0.000000e+00},
  Vector3:new{x= 3.000000e-02, y=0.000000e+00},
  Vector3:new{x= 1.200000e-01, y=0.000000e+00},
  Vector3:new{x= 1.600000e-01, y=0.000000e+00},
  Vector3:new{x= 2.600000e-01, y=0.000000e+00},
  Vector3:new{x= 4.000000e-01, y=0.000000e+00},
  Vector3:new{x= 5.700000e-01, y=0.000000e+00},
  Vector3:new{x= 7.700000e-01, y=0.000000e+00},
  Vector3:new{x= 1.000000e+00, y=0.000000e+00},
}

r = ctrl_pts[1].y
bd = Bezier:new{points=ctrl_pts}
a = Vector3:new{x=ctrl_pts[1].x,  y=0.0}
b = Vector3:new{x=ctrl_pts[1].x,  y=ctrl_pts[1].y}
c = Vector3:new{x=ctrl_pts[npoints].x, y=0.0}
d = Vector3:new{x=ctrl_pts[npoints].x, y=ctrl_pts[npoints].y}
ab = Line:new{p0=a, p1=b}
--ac = Line:new{p0=a, p1=c}
ac = Bezier:new{points=ctrl_pts2}
cd = Line:new{p0=c, p1=d}
aa = Vector3:new{x=a.x - r, y=0.0}
bb = Vector3:new{x=b.x - r, y=b.y}
aaa = Line:new{p0=aa, p1=a}
aabb= Line:new{p0=aa, p1=bb}
bbb = Line:new{p0=bb, p1=b}
L_noz = ctrl_pts[npoints].x - ctrl_pts[1].x
expansion_ratio = c.y/a.y

throat_patch = CoonsPatch:new{north=bbb, east=ab, south=aaa, west=aabb}
expand_patch = CoonsPatch:new{north=bd, east=cd, south=ac, west=ab}

-- discretisation
nni = 460
nnj = 64
nni_t = 32

grid = StructuredGrid:new{psurface=expand_patch, niv=nni, njv=nnj,
    cfList = {north=QuadraticFunction:new{ratio=2.0}, south=QuadraticFunction:new{ratio=2.0},
              west=GeometricFunction:new{a=0.00005, r=1.3, N=nnj, reverse=true},
              east=GeometricFunction:new{a=0.00005, r=1.3, N=nnj, reverse=true}
          }
}
registerFluidGridArray{
   grid=grid,
   nib=ncores-1, njb=1,
   fsTag='initial',
   bcTags={north='wall', east='outflow'}
}

throatgrid = StructuredGrid:new{psurface=throat_patch, niv=nni_t, njv=nnj,
    cfList = {
              west=GeometricFunction:new{a=0.00005, r=1.3, N=nnj, reverse=true},
              east=GeometricFunction:new{a=0.00005, r=1.3, N=nnj, reverse=true}
          }
}

registerFluidGridArray{
   grid=throatgrid,
   nib=1, njb=1,
   fsTag='initial',
   bcTags={west='inflow', north='wall'}
}
identifyGridConnections()


flowDict = {
   initial=initial,
}

bcDict = {
    inflow=InFlowBC_Supersonic:new{flowState=inflow},
    wall=WallBC_NoSlip_FixedT:new{Twall=300.0, group="wall"},
    outflow=OutFlowBC_Simple:new{}
}

makeFluidBlocks(bcDict, flowDict)

turbZone_botLeftCorner = Vector3:new{x=0.10, y=y0}
turbZone_topRightCorner = Vector3:new{x=1.2, y=0.2}
TurbulentZone:new{p0=turbZone_botLeftCorner, p1=turbZone_topRightCorner}

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 4000,
   stop_on_relative_residual = 1.0e-8,
   number_of_phases = 2,
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = true,
   max_linear_solver_iterations = 40,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,
   write_loads = true,
   total_snapshots = 2,
   steps_between_snapshots = 100,
   steps_between_diagnostics = 1,
   steps_between_loads_update = 1000000,
   max_steps_in_initial_phases = { 2000, },
   phase_changes_at_relative_residual = { 1.0e-06, },
   allowable_relative_mass_change = 0.2,
   min_relaxation_factor_for_update = 0.0001,
   min_relaxation_factor_for_cfl_growth = 0.5,
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   frozen_preconditioner = true,
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 10,
   linear_solve_tolerance = 0.01,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.1,
   start_cfl = 0.5,
   max_cfl = 1.0e4,
   auto_cfl_exponent = 0.7,
   limit_on_cfl_decrease_ratio = 1.0,
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 5.0,
   max_cfl = 1.0e4,
}
