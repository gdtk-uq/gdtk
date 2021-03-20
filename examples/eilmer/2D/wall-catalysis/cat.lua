-- Blunt wedge with catalytic wall boundary. Geometry from:
--
-- Hicham Alkandry, Kyle Hanquist, and Iain D. Boyd. 
-- "Conceptual Analysis of Electron Transpiration Cooling 
-- for the Leading Edges of Hypersonic Vehicles", 
-- 11th AIAA/ASME Joint Thermophysics and Heat Transfer Conference, 
-- AIAA AVIATION Forum, (AIAA 2014-2674) 
-- https://doi.org/10.2514/6.2014-2674 
--
-- @author: Nick N. Gibbons
-- 2021-01-21

job_title = "Blunt Wedge with catalytic wall"
print(config.title)
config.dimensions = 2
config.axisymmetric = false
config.viscous = true
config.viscous_signal_factor = 1.0
config.spatial_deriv_locn = "cells"
config.spatial_deriv_calc = "least_squares"
config.mass_diffusion_model = "ficks_first_law"
config.diffusion_coefficient_type = "binary_diffusion"

N_flows = 20
N_solutions = 20
flow_time = 0.02 / 6000 -- 2 radaii / flow velocity
config.max_time = N_flows*flow_time
config.dt_plot = config.max_time/N_solutions
config.dt_init = 1.0e-13
config.cfl_value = 0.50
config.max_step = 100000000
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "rkl1"
config.interpolation_delay = 0.1*flow_time -- Startup phase can be crashy without this
config.print_count=100

config.write_loads = true
config.dt_loads = config.max_time/N_solutions
config.boundary_groups_for_loads = "wall"

nsp, nmodes, gmodel = setGasModel('air-5sp-2T.lua')
config.reacting = true
config.reactions_file = 'air-5sp-2T-chem.lua'
config.energy_exchange_file = 'air-energy-exchange.lua'

T_inf = 238.0  -- K
rho_inf = 2.3e-4  -- kg/m^3
mass_fraction = {N2=0.767, O2=0.233}
u_inf = 6000.0  -- m/s

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

initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0.0, vely=0.0, massf=mass_fraction, T_modes={T_inf}}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction, T_modes={T_inf}}

-- Grid Geometry Specs  --            bo   ___---* co
ri = 0.010              --        _-- *---/      |    -- Notes: - ai is always at the origin
thetai = math.rad(5.0)  --       /               |              - ci and co are at x=L
ro = 0.0215             --      /    _-*---------* ci           - oo is downstream of oi by diffo
diffo = 0.0075          --     /    /  bi           
thetao = math.rad(45.0) --     |    |
L = 0.018               --  ao *----*  *  *

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
ao = aoco(0.0)
co = aoco(1.0)

aoai = Line:new{p0=ao, p1=ai}
coci = Line:new{p0=co, p1=ci}

niv =40+1; njv =40+1;
surface = CoonsPatch:new{north=coci, south=aoai, east=aici, west=aoco}
cluster_south = GeometricFunction:new{a=0.005, r = 1.2, N=niv, reverse=true}
cluster_north = GeometricFunction:new{a=0.005, r = 1.2, N=niv, reverse=true}
clusterlist = {north=cluster_north, south=cluster_south, east=none, west=none}

grid = StructuredGrid:new{psurface=surface, niv=niv, njv=njv, cfList=clusterlist}

-- Three different kinds of catalytic behaviour are available,
-- Uncomment one of the "east" lines to try them out.
blocks = FluidBlockArray{grid=grid,
                        initialState=initial,
                        bcList={south=WallBC_WithSlip:new{},
                                --east=WallBC_NoSlip_FixedT:new{catalytic_type="equilibrium", Twall=2200.0, group="wall"},
                                east=WallBC_NoSlip_UserDefinedT:new{catalytic_type="equilibrium", Twall="udf-twall.lua", group="wall"},
                                --east=WallBC_NoSlip_FixedT:new{catalytic_type="none", Twall=2200.0, group="wall"},
                                --east=WallBC_NoSlip_FixedT:new{catalytic_type="fixed_composition", wall_massf_composition=mass_fraction, Twall=2200.0, group="wall"},
                                north=OutFlowBC_Simple:new{},
                                west=InFlowBC_Supersonic:new{flowState=inflow}},
                        nib=4, njb=4}


identifyBlockConnections()
