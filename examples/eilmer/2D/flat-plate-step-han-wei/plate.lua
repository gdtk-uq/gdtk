-- plate.lua
--
-- Initial code:
-- Han Wei 19-Feb-2017
--
-- History:
-- The grid is revised, mainly to increase the y coordinate
-- of the downstream end of the upper boundary as it's found
-- to be too small for the fully heated case (19-Oct-2017).
-- The history locations are updated with the actual sensor locations
-- in the final model design, heat flux data is output for comparisons,
-- and the simulation time is extended to 400us.
--
-- Eilmer4 example:
-- PJ 2018-03-17 ported from Han's Python input script
-- and brought over 5-species air chemistry from Rowan's sphere-lobb case.

config.title = "X2 cold-1 flow condition over a flat plate with "..
   "a 8mm high step with the model unheated"
print(config.title)
config.dimensions = 2
nsp, nmodes, gm = setGasModel('air-5sp.lua')
config.reacting = true
config.reactions_file = 'air-5sp-6r.lua'

-- Flow conditions
p_inf = 349.0  -- Pa
u_inf = 3890.0 -- m/s
T_inf = 961.0  -- K
p_ini = 40.0   -- Pa
T_ini = 298.15 -- K
T_wall = 298.15 -- K
T_heating = T_wall -- Temperature of the heated section, K
T_pl = T_wall -- Temperature of the rest of the plate, K

-- mole fraction for the inflow from CEA calculation
molefif= {N2=7.8999e-1, NO=2.0045e-5, O2=2.0999e-1, N=0.0, O=3.555e-10}
massfif = gm:molef2massf(molefif)

initial = FlowState:new{p=p_ini, T=T_ini, massf={N2=0.767, O2=0.233}}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, massf=massfif}

-- Making the grid
-- Geometry
lle = 0.018 -- length of the leading edge, m
lplate = 0.100 -- length of the flat plate, m
lheated = 0.045 -- length of the heated section, m
lcold = 0.037 -- length of the unheated section, m
lsp = 0.018 -- length of the estimated separation region, m
hstep = 0.008 -- height of the step, m
lstep = 0.010 -- length of the step, m
dedge = 3.0*hstep+0.003
-- distance between the edge of the computational domain
-- and the surface of the plate
d0 = hstep
-- y distance between the usptream end of the computational domain
-- and the surface of the plate

yv = function(x, x1, y1, x2, y2)
   x1 = x1 or 0.0
   y1 = y1 or d0
   x2 = x2 or (lle+lheated+lcold+lstep)
   y2 = y2 or dedge
   return ((y2-y1)/(x2-x1))*(x-x1)+y1
end

yv0 = function(x, x1, y1, x2, y2)
   x1 = x1 or 0.0
   y1 = y1 or (0.4*hstep)
   x2 = x2 or (lle+lheated+lcold)
   y2 = y2 or hstep
   return ((y2-y1)/(x2-x1))*(x-x1)+y1
end

pnts = {
   a = Vector3:new{x=0.0, y=0.0},
   b = Vector3:new{x=lle, y=0.0},
   c = Vector3:new{x=lle+lheated, y=0.0},
   d = Vector3:new{x=lle+lheated+lcold-lsp, y=0.0},
   e = Vector3:new{x=lle+lheated+lcold, y=0.0},
   f = Vector3:new{x=0.0, y=0.4*hstep},
   g = Vector3:new{x=lle, y=yv0(lle)},
   h = Vector3:new{x=lle+lheated, y=yv0(lle+lheated)},
   i = Vector3:new{x=lle+lheated+lcold-lsp, y=yv0(lle+lheated+lcold-lsp)},
   j = Vector3:new{x=lle+lheated+lcold, y=hstep},
   k = Vector3:new{x=lle+lheated+lcold+lstep, y=hstep},
   l = Vector3:new{x=0.0, y=d0},
   m = Vector3:new{x=lle, y=yv(lle)},
   n = Vector3:new{x=lle+lheated, y=yv(lle+lheated)},
   o = Vector3:new{x=lle+lheated+lcold-lsp, y=yv(lle+lheated+lcold-lsp)},
   p = Vector3:new{x=lle+lheated+lcold, y=yv(lle+lheated+lcold)},
   q = Vector3:new{x=lle+lheated+lcold+lstep, y=dedge},
   r = Vector3:new{x=lle+lheated+lcold+lstep, y=0.0}
}


blkLabels = {[0]="Surface-1", [1]="Surface-2", [2]="Surface-3", [3]="Surface-4",
   [4]="Edge-1", [5]="Edge-2", [6]="Edge-3", [7]="Edge-4", [8]="Edge-5"}

bcLists = {}
bcLists[0] = {north=nil, east=nil,
              south=WallBC_NoSlip_FixedT:new{Twall=T_wall, label="loads"}, 
              west=InFlowBC_Supersonic:new{flowState=inflow}}
bcLists[1] = {north=nil, east=nil, west=nil,
              south=WallBC_NoSlip_FixedT:new{Twall=T_heating, label="loads"}}
bcLists[2] = {north=nil, east=nil, west=nil,
              south=WallBC_NoSlip_FixedT:new{Twall=T_heating, label="loads"}}
bcLists[3] = {north=nil, west=nil,
              east=WallBC_NoSlip_FixedT:new{Twall=T_wall, label="loads"},
              south=WallBC_NoSlip_FixedT:new{Twall=T_heating, label="loads"}}
bcLists[4] = {north=InFlowBC_Supersonic:new{flowState=inflow},
              east=nil, south=nil,
              west=InFlowBC_Supersonic:new{flowState=inflow}}
bcLists[5] = {north=InFlowBC_Supersonic:new{flowState=inflow},
              east=nil, south=nil, west=nil}
bcLists[6] = {north=InFlowBC_Supersonic:new{flowState=inflow},
              east=nil, south=nil, west=nil}
bcLists[7] = {north=InFlowBC_Supersonic:new{flowState=inflow},
              east=nil, south=nil, west=nil}
bcLists[8] = {north=InFlowBC_Supersonic:new{flowState=inflow},
              east=OutFlowBC_Simple:new{},
              south=nil, west=nil}

-- p01--p11
--  |    |
-- p00--p10
quads = {}
-- Start with lower layer of quads, along plate surface.
quads[0] = CoonsPatch:new{p01=pnts.f, p11=pnts.g, p10=pnts.b, p00=pnts.a}
quads[1] = CoonsPatch:new{p01=pnts.g, p11=pnts.h, p10=pnts.c, p00=pnts.b}
quads[2] = CoonsPatch:new{p01=pnts.h, p11=pnts.i, p10=pnts.d, p00=pnts.c}
quads[3] = CoonsPatch:new{p01=pnts.i, p11=pnts.j, p10=pnts.e, p00=pnts.d}
-- Upper layer of quads.
quads[4] = CoonsPatch:new{p01=pnts.l, p11=pnts.m, p10=pnts.g, p00=pnts.f}
quads[5] = CoonsPatch:new{p01=pnts.m, p11=pnts.n, p10=pnts.h, p00=pnts.g}
quads[6] = CoonsPatch:new{p01=pnts.n, p11=pnts.o, p10=pnts.i, p00=pnts.h}
quads[7] = CoonsPatch:new{p01=pnts.o, p11=pnts.p, p10=pnts.j, p00=pnts.i}
quads[8] = CoonsPatch:new{p01=pnts.p, p11=pnts.q, p10=pnts.k, p00=pnts.j}

-- Set calculation size
fullRes = false
-- Basic discretization for a block.
if fullRes then
   nx = 25; ny = 25 -- Han's full resolution
else
   nx = 15; ny = 15 -- cut-down exercise
end
-- Numbers of blocks within block arrays.
nib = {}; njb = {}
if fullRes then -- Han's full resolution
   nib[0] =  7; njb[0] = 6
   nib[1] = 14; njb[1] = 6
   nib[2] =  7; njb[2] = 6
   nib[3] =  6; njb[3] = 6
   nib[4] =  7; njb[4] = 4
   nib[5] = 14; njb[5] = 4
   nib[6] =  7; njb[6] = 4
   nib[7] =  6; njb[7] = 4
   nib[8] =  2; njb[8] = 4
else -- cut-down exercise
   nib[0] =  4; njb[0] = 3
   nib[1] =  7; njb[1] = 3
   nib[2] =  4; njb[2] = 3
   nib[3] =  3; njb[3] = 3
   nib[4] =  4; njb[4] = 2
   nib[5] =  7; njb[5] = 2
   nib[6] =  4; njb[6] = 2
   nib[7] =  3; njb[7] = 2
   nib[8] =  1; njb[8] = 2
end

-- Grid clustering
clb = RobertsFunction:new{end0=true,end1=false,beta=1.3}
clb_med = RobertsFunction:new{end0=true,end1=false,beta=1.5}
clb_lar = RobertsFunction:new{end0=true,end1=false,beta=1.2}
crb = RobertsFunction:new{end0=false,end1=true,beta=1.3}
crb_med = RobertsFunction:new{end0=false,end1=true,beta=1.6}
cfLists = {}
cfLists[0] = {north=clb, east=clb_med, south=clb, west=clb_med}
cfLists[1] = {north=crb_med, east=clb_med, south=crb_med, west=clb_med}
cfLists[2] = {north=nil, east=clb_med, south=nil, west=clb_med}
cfLists[3] = {north=nil, east=clb_med, south=nil, west=clb_med}
cfLists[4] = {north=clb, east=clb_med, south=clb, west=clb_med}
cfLists[5] = {north=crb_med, east=clb_med, south=crb_med, west=clb_med}
cfLists[6] = {north=nil, east=clb_med, south=nil, west=clb_med}
cfLists[7] = {north=nil, east=clb_med, south=nil, west=clb_med}
cfLists[8] = {north=clb_lar, east=clb_med, south=clb_lar, west=clb_med}

grids = {}
for i=0,8 do
   grids[i] = StructuredGrid:new{psurface=quads[i], cfList=cfLists[i],
                                 niv=nx*nib[i]+1, njv=ny*njb[i]+1}
end

blks = {}
for i=0,8 do
   blks[i] = FBArray:new{grid=grids[i], nib=nib[i], njb=njb[i], 
                             initialState=initial, bcList=bcLists[i],
                             label=blkLabels[i]}
end
identifyBlockConnections()

-- Do a little more setting of global data.
config.gasdynamic_update_scheme = "predictor-corrector"
config.viscous = true
config.flux_calculator = "adaptive"
config.max_time = 4.0e-4  -- seconds
config.max_step = 2500000
config.dt_init = 1.0e-10
-- config.dt_max = 1.0e-9
config.cfl_value = 1.0
config.dt_plot = config.max_time / 200.0
config.dt_history = 1.0e-6

-- PJ wants to run the calculation on his 8-core workstation.
mpiDistributeBlocks{ntasks=8, dist="load-balance"}

dsen = 0.0045 -- distance between each sensor, m
-- Kulite pressure transducer or thin film heat transfer gauge
setHistoryPoint{x=lplate-11.5*dsen, y=0.0} -- Sensor-01
setHistoryPoint{x=lplate-10.5*dsen, y=0.0} -- Sensor-02
setHistoryPoint{x=lplate-9.5*dsen, y=0.0} -- Sensor-03
setHistoryPoint{x=lplate-8.5*dsen, y=0.0} -- Sensor-04
setHistoryPoint{x=lplate-7.5*dsen, y=0.0} -- Sensor-05
setHistoryPoint{x=lplate-6.5*dsen, y=0.0} -- Sensor-06
setHistoryPoint{x=lplate-5.5*dsen, y=0.0} -- Sensor-07
setHistoryPoint{x=lplate-4.5*dsen, y=0.0} -- Sensor-08
setHistoryPoint{x=lplate-3.5*dsen, y=0.0} -- Sensor-09
setHistoryPoint{x=lplate-2.5*dsen, y=0.0} -- Sensor-10
setHistoryPoint{x=lplate-1.5*dsen, y=0.0} -- Sensor-11
setHistoryPoint{x=lplate-0.5*dsen, y=0.0} -- Sensor-12
-- These two sensor locations fall on the step
setHistoryPoint{x=lplate+0.5*dsen, y=0.0} -- Sensor-13
setHistoryPoint{x=lplate+1.5*dsen, y=0.0} -- Sensor-14

dofile("sketch-domain.lua")
