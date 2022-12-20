-- Simulation of ablation experiments, loosely based on:
-- "Nonequilibrium Finite-Rate Carbon Ablation Model for Earth Reentry Flows"
-- Christopher R. Alba and Robert B. Greendyke, JSR Vol. 53, No. 3, May–June 2016
-- 
-- @author: Nick N. Gibbons 
-- 2021-06-25

job_title = "Ablation Stuff"
config.dimensions = 2
config.axisymmetric = true 
config.flux_calculator = "hanel"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "shock_fitting"

nsp, nmodes, gmodel = setGasModel('gm-air11-2T.lua')
config.reacting = true
config.reactions_file = 'rr-kim-air11-2T.lua'
config.energy_exchange_file = 'ee-kim-air11-2T.lua'

u_inf = 7000.0
rho_inf = 0.018011932984430126
T_inf = 226.65
mass_fraction = {N2=0.767, O2=0.233}

-- Compute full inflow state from T,rho, and u
Q = GasState:new{gmodel}
Q.T = T_inf;
Q.rho = rho_inf;
Q.massf = mass_fraction
Q.T_modes = {T_inf}
gmodel:updateThermoFromRHOT(Q);
p_inf = Q.p

initial = FlowState:new{p=p_inf/2, T=T_inf, velx=u_inf/2, vely=0.0, massf=mass_fraction, T_modes={T_inf}}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction, T_modes={T_inf}}

-- Grid Geometry Specs  --            bo   ___---* co
ri = 0.010              --        _-- *---/      |    -- Notes: - ai is always at the origin
thetai = math.rad(8.0)  --       /               |              - ci and co are at x=L
ro = 0.015              --      /    _-*---------* ci           - oo is downstream of oi by diffo
diffo = 0.003           --     /    /  bi           
thetao = math.rad(30.0) --     |    |
L = 0.0202+ri           --  ao *----*ai *oi *oo

oi = Vector3:new{x=ri, y=0.0}
ai = oi + ri*Vector3:new{x=-1.0, y=0.0}
bi = oi + ri*Vector3:new{x=-math.sin(thetai), y=math.cos(thetai)}
dxci = L - bi.x
dyci = dxci*math.tan(thetai)
ci = Vector3:new{x=bi.x + dxci, y=bi.y+dyci}
aibi = Arc:new{p0=ai, p1=bi, centre=oi}
bici = Line:new{p0=bi, p1=ci}
aici = Polyline:new{segments={aibi, bici}}

oo = Vector3:new{x=ri+diffo, y=0.0}
ao = oo + ro*Vector3:new{x=-1.0, y=0.0}
bo = oo + ro*Vector3:new{x=-math.sin(thetao), y=math.cos(thetao)}
dxco = L - bo.x
dyco = dxco*math.tan(thetao)
co = Vector3:new{x=bo.x + dxco, y=bo.y+dyco}
aobo = Arc:new{p0=ao, p1=bo, centre=oo}
boco = Line:new{p0=bo, p1=co}
aoco = Polyline:new{segments={aobo, boco}}

aoai = Line:new{p0=ao, p1=ai}
coci = Line:new{p0=co, p1=ci}

niv =32+1; njv =86+1;
surface = CoonsPatch:new{north=coci, south=aoai, east=aici, west=aoco}

cluster_south = GaussianFunction:new{m=0.01, s=0.15,ratio=0.1}
cluster_north = GaussianFunction:new{m=0.01, s=0.15,ratio=0.1}
cluster_east = GaussianFunction:new{m=0.01, s=0.8, ratio=0.1}
cluster_west = GaussianFunction:new{m=0.01, s=0.8, ratio=0.1}
clusterlist = {north=cluster_north, south=cluster_south, east=cluster_east, west=cluster_west}

grid = StructuredGrid:new{psurface=surface, niv=niv, njv=njv, cfList=clusterlist}

blocks = FBArray:new{grid=grid,
                        initialState=initial,
                        bcList={
                                north=OutFlowBC_Simple:new{},
                                west=InFlowBC_ShockFitting:new{flowState=inflow}},
                        nib=1, njb=8}



identifyBlockConnections()

N_flows = 14
N_solutions = 2
flow_time = ri/u_inf -- 2 radaii / flow velocity
config.max_time = N_flows*flow_time
config.dt_plot = config.max_time/N_solutions
config.dt_init = 5.0e-14
config.cfl_value = 0.5
config.max_step = 100000000
config.print_count=100
config.control_count=100

config.shock_fitting_delay = 2.0*flow_time
config.interpolation_delay = 6.0*flow_time
