-- Simulation of Fire II Flight Experiments
-- From doi.org/10.2514/6.2007-605
-- @author: Nick N. Gibbons 
-- 2020-11-30

job_title = "Fire II Inviscid Flow"
config.dimensions = 2
config.axisymmetric = true
config.flux_calculator = "hanel"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "shock_fitting"

nsp, nmodes, gmodel = setGasModel('gm-air11-2T.lua')
config.reacting = true
config.reactions_file = 'rr-park-air11-2T.lua'
config.energy_exchange_file = 'ee-park-air11-2T.lua'

-- 1636 Second Trajectory Point, at 71.04 km Altitude
T_inf = 210      -- K
T_wall= 810
rho_inf = 8.57e-5  -- kg/m^3
mass_fraction = {N2=0.767, O2=0.233}
u_inf = 11.31e3    -- m/s

-- Compute full inflow state from T,rho, and u
Q = GasState:new{gmodel}
Q.T = T_inf;
Q.rho = rho_inf;
Q.massf = mass_fraction
Q.T_modes = {T_inf}
gmodel:updateThermoFromRHOT(Q);
gmodel:updateSoundSpeed(Q)
p_inf = Q.p
M_inf = u_inf/Q.a
gmodel:updateTransCoeffs(Q);

initial = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction, T_modes={T_inf}}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction, T_modes={T_inf}}

-- Grid Geometry Specs  
Ri = 0.9347             
ri = 0.0102
A  = 0.3358
thetai = math.asin((A-ri)/(Ri-ri))
L = 0.05
diffo = 0.07
Ro = Ri + diffo
ro = ri + diffo
psio = math.rad(5.0)

oi = Vector3:new{x=Ri, y=0.0}
ai = oi + Ri*Vector3:new{x=-1.0, y=0.0}
bi = oi + Ri*Vector3:new{x=-math.cos(thetai), y=math.sin(thetai)}
pi = oi + (Ri-ri)*Vector3:new{x=-math.cos(thetai), y=math.sin(thetai)}
ci = pi + ri*Vector3:new{x=0.0, y=1.0}
di = ci + L*Vector3:new{x=1.0, y=0.0}

aibi = Arc:new{p0=ai, p1=bi, centre=oi}
bici = Arc:new{p0=bi, p1=ci, centre=pi}
cidi = Line:new{p0=ci, p1=di}
aidi = Polyline:new{segments={aibi, bici, cidi}}

thetao = 1.5*thetai
oo = Vector3:new{x=Ri, y=0.0}
ao = oo + Ro*Vector3:new{x=-1.0, y=0.0}
do_= oo + Ro*Vector3:new{x=-math.cos(thetao), y=math.sin(thetao)}

aodo = Arc:new{p0=ao, p1=do_, centre=oo}
aoai = Line:new{p0=ao, p1=ai}
dodi = Line:new{p0=do_, p1=di}

--niv =64+1; njv =256+1;
niv =32+1; njv =128+1;
surface = CoonsPatch:new{north=dodi, south=aoai, east=aidi, west=aodo}

cluster_east = GaussianFunction:new{m=0.001, s=0.7, ratio=0.5}
cluster_west = GaussianFunction:new{m=0.001, s=0.5, ratio=0.5}
clusterlist = {north=none, south=none, east=cluster_east, west=cluster_west}

grid = StructuredGrid:new{psurface=surface, niv=niv, njv=njv, cfList=clusterlist}

blocks = FBArray:new{grid=grid,
                        initialState=initial,
                        bcList={
                                north=OutFlowBC_Simple:new{},
                                west=InFlowBC_ShockFitting:new{flowState=inflow}},
                        nib=2, njb=8}




identifyBlockConnections()

N_flows = 12.0
N_solutions = 24.0
flow_time = Ri/u_inf
config.max_time = N_flows*flow_time
config.dt_plot = config.max_time/N_solutions
config.dt_init = 1.0e-16
config.cfl_value = 0.2
config.max_step = 2147483647
config.print_count=100
config.shock_fitting_delay = 1.0*flow_time
config.shock_fitting_scale_factor = 0.1
config.flowstate_limits_max_temp = 100000
config.interpolation_delay = 4.0*flow_time
config.allow_reconstruction_for_species = false
