-- sphere.lua
--
-- X3 Argon Experimental test flow over spherical model.
-- Daniel Smith, November 2019
-- PJ, 2019-12-31, adjustments for Eilmer4
--     2020-01-17, use Daniel's 2T gas model

config.title = "X3_Finite_Rate_Argon"
config.axisymmetric = true
config.flux_calculator = "hanel"
config.extrema_clipping = false
config.interpolation_order=1

nsp, nmodes, gm = setGasModel("argon-gas-2T.lua")
config.reacting = true
config.reactions_file = "argon-gas-2T.lua"

-- Define flow conditions
p_inf = 342.0 -- Pa
T_inf = 280.0 -- degrees K
v_inf = 9000.0 -- flow speed, m/s
mole_inf = {['Ar']=0.98, ['Ar+']=0.01, ['e-']=0.01}
mf_inf = gm:molef2massf(mole_inf)

inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf}, velx=v_inf, massf=mf_inf}
initial = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf}, velx=v_inf, massf=mf_inf}
--initial = FlowSolution:new{jobName="sphere", dir="../noreactions", tindx=31, nBlocks=2}

--geometry
factor = 2.0
nnEW = 30*factor                        
nnNS0 = 80*factor
nnNS1 = 40*factor

R = 0.0325      -- m
delta = 0.0040  --distance to front of mesh
-- (original was 0.0065 but 0.0080 needed for lower speed sims)
h2 = 0.057      --height of mesh at exit
h1 = 0.027      --height of mesh at shoulder
d = 0.05        --length of straight section

pnts = {}
pnts.Cntr = {x=0, y= 0}
pnts.A = {x=-R, y= 0.0}
pnts.B = {x=0.0, y= R}
pnts.C = {x=d, y= R}
pnts.D = {x=d, y= R+h2}
pnts.E = {x=0.0, y= R+h1}
pnts.F = {x=-R-delta, y= 0.0}
pnts.Bez1 = {x=-R-delta, y=pnts.E.y-(pnts.E.x-(-R-delta))*(h2-h1)/d}

east0  = Arc:new{p0=pnts.A, p1=pnts.B, centre=pnts.Cntr}
north0 = Line:new{p0=pnts.E, p1=pnts.B}
west0  = Bezier:new{points={pnts.F, pnts.Bez1, pnts.E}}
south0 = Line:new{p0=pnts.F, p1=pnts.A}

quad0 = CoonsPatch:new{north=north0, east=east0, south=south0, west=west0}
quad1 = CoonsPatch:new{p00=pnts.E, p10=pnts.B, p11=pnts.C, p01=pnts.D}

cf = {east = GaussianFunction:new{m=0.01, s=0.8, ratio=0.1}, 
      west = GaussianFunction:new{m=0.01, s=0.8, ratio=0.1}} 

grid0 = StructuredGrid:new{psurface=quad0, niv=nnEW, njv=nnNS0, cfList=cf}
grid1 = StructuredGrid:new{psurface=quad1, niv=nnEW, njv=nnNS1}

blk0 = FBArray:new{grid=grid0, initialState=initial, nib=2, njb=4, bcList={west=InFlowBC_Supersonic:new{flowState=inflow}}}
blk1 = FBArray:new{grid=grid1, initialState=initial, nib=2, njb=2, bcList={west=InFlowBC_Supersonic:new{flowState=inflow},north=OutFlowBC_Simple:new{}}}

identifyBlockConnections()

--Config options
SteadyStateSolver{
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 25,
   
   use_complex_matvec_eval = true,
   number_total_steps = 2000,
   stop_on_relative_global_residual = 1e-9,
   use_physicality_check = true,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 80,
   max_restarts = 0,
   residual_based_cfl_scheduling = true,
   cfl_max = 1e4,

   -- Settings for start-up phase
   number_start_up_steps = 99,
   cfl0 = 0.1,
   eta0 = 0.001,
   tau0 = 1e-2,
   sigma0 = 1.0e-30,
   p0 = 0.7,

   -- Settings for inexact Newton phase
   cfl1 = 1.0,
   eta1 = 0.01,
   tau1 = 5e-1,
   sigma1 = 1.0e-30,
   eta_strategy = "constant",
   p1 = 0.7,

   -- Settings control write-out
   snapshots_count = 100,
   number_total_snapshots = 100,
   write_diagnostics_count = 1,
   write_loads_count = 100,
}

