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
config.viscous = true
config.mass_diffusion_model = "ficks_first_law"
config.diffusion_coefficient_type = "binary_diffusion"
config.enforce_species_density_positivity = true
config.interpolation_order=1

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

aoco = Spline2:new{filename="../shockfit/inflated_shock_shape.dat"}
ao = aoco(0.0)
co = aoco(1.0)

aoai = Line:new{p0=ao, p1=ai}
coci = Line:new{p0=co, p1=ci}

--niv =64+1; njv =120+1;
niv =48+1; njv =86+1;
surface = CoonsPatch:new{north=coci, south=aoai, east=aici, west=aoco}

cluster_south = GaussGeomHybridFunction:new{A=0.001, R=1.2, N=niv, m=0.65, s=0.20,ratio=0.2, reverse=true}
cluster_north = GaussGeomHybridFunction:new{A=0.001, R=1.2, N=niv, m=0.65, s=0.20,ratio=0.2, reverse=true}
cluster_east = GaussianFunction:new{m=0.01, s=0.8, ratio=0.1}
cluster_west = GaussianFunction:new{m=0.01, s=0.8, ratio=0.1}
clusterlist = {north=cluster_north, south=cluster_south, east=cluster_east, west=cluster_west}

grid = StructuredGrid:new{psurface=surface, niv=niv, njv=njv, cfList=clusterlist}

blocks = FBArray:new{grid=grid,
                        initialState=initial,
                        bcList={
                                north=OutFlowBC_Simple:new{},
                                west=InFlowBC_Supersonic:new{flowState=inflow},
                                east=WallBC_NoSlip_FixedT:new{Twall=2000.0, group="wall",
                                                              user_post_diff_flux="ablation-udf.lua",
                                   }
                                },
                        nib=1, njb=8}



identifyBlockConnections()

config.write_loads = true
config.boundary_groups_for_loads = "wall"

flowtime = L/u_inf
config.gasdynamic_update_scheme = "backward_euler"
config.cfl_schedule = {{flowtime*0.1,1.0},
                       {flowtime*0.5,1.0},
                       {flowtime*1.0,2.0},
                       {flowtime*2.0,20.0}}
config.max_time = flowtime*3.0  -- About 4 flow lengths (L=1.4m, 1 flow length ~ 1.96 ms)
config.dt_plot = config.max_time/20.0 
config.max_step = 300000

