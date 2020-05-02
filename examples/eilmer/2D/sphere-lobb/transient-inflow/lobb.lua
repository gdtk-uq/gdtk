-- Authors: Rowan G. and Peter J.
-- Date: 2016-01-19, simplified 2020-05-02
--
-- Description:
--   A nylon sphere, 0.5-inch diameter, fired into a ballistic test range.
--
-- Reference:
--  Lobb, R.K. (1964)
--  Experimental measurement of shock detachment distance on spheres
--  fired in air at hypervelocities.
--  pp 519--527 in
--  The High Temperature Aspects of Hypersonic Flow
--  edited by Nelson, W.C.
--  Pergamon Press, Oxford, 1964

config.title = "Sphere fired into air."
print(config.title)
Db = 0.5 * 0.0254 -- diameter (in m) of ball bearing
Rc = Db/2

setGasModel('air-5sp.lua')
p_inf = 666.0 -- Pa (=5 Torr)
T_inf = 293.0 -- K
vx_inf = 4825.0 -- m/s
initial = FlowState:new{p=p_inf, T=T_inf, velx=1000.0,
                        massf={N2=0.78,O2=0.22}}
M_inf = vx_inf/initial.a -- used for Billig correlation
body_flow_time = Db/vx_inf

config.reacting = true
config.reactions_file = 'air-5sp-6r.lua'
config.axisymmetric = true
config.max_time = 12*body_flow_time
print("max_time=", config.max_time)
config.max_step = 80000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/10

require "billig_patch"
bp = billig_patch.make_patch{Minf=M_inf, R=Rc, scale=0.9}
cf_circum = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grd = StructuredGrid:new{psurface=bp.patch, niv=61, njv=61,
                         cfList={west=cf_circum, east=cf_circum}}
blk = FluidBlockArray{
   grid=grd, initialState=initial, label='blk',
   bcList={west=InFlowBC_Transient:new{filename='inflow.data'},
           north=OutFlowBC_Simple:new{}},
   nib=1, njb=4}

setHistoryPoint{x=-Rc, y=0.0} -- stagnation point
setHistoryPoint{x=0.0, y=Rc}  -- 90 degree point on sphere
config.dt_history = 1.0e-7
