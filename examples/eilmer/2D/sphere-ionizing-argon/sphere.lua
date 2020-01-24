-- sphere.lua
--
-- X3 Argon Experimental test flow over spherical model.
-- Daniel Smith, November 2019
-- PJ, 2019-12-31, adjustments for Eilmer4
--     2020-01-17, use Daniel's 2T gas model

config.title = "X3_Finite_Rate_Argon"
config.axisymmetric = true

nsp, nmodes, gm = setGasModel("argon-gas-2T.lua")
config.reacting = true
config.reactions_file = "argon-gas-2T.lua"

-- Define flow conditions
p_inf = 342.0 -- Pa
T_inf = 280.0 -- degrees K
v_inf = 5728.0 -- flow speed, m/s
mf_inf = {['Ar']=1.0, ['Ar+']=0.0, ['e-']=0.0}

inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf},
                       velx=v_inf, massf=mf_inf}
initial = FlowState:new{p=0.3*p_inf, T=T_inf, T_modes={T_inf},
                        velx=0.6*v_inf, massf=mf_inf}

--geometry
nnEW = 30                         
nnNS0 = 60
nnNS1 = 50

R = 0.0325      -- m
delta = 0.0065  --distance to front of mesh
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

grid0 = StructuredGrid:new{psurface=quad0, niv=nnEW, njv=nnNS0}
grid1 = StructuredGrid:new{psurface=quad1, niv=nnEW, njv=nnNS1}

blk0 = FluidBlock:new{grid=grid0, initialState=initial}
blk1 = FluidBlock:new{grid=grid1, initialState=initial}

identifyBlockConnections()

blk0.bcList[west] = InFlowBC_Supersonic:new{flowState=inflow}
blk1.bcList[west] = InFlowBC_Supersonic:new{flowState=inflow}
blk1.bcList[north] = OutFlowBC_Simple:new{}

--Config options
config.sticky_electrons = true -- just wanted to try it
config.cfl_value = 0.05 -- to get better chemistry-gas-dynamics coupling

config.max_time = 100.0e-6
config.max_step = 200000
config.dt_plot = config.max_time/50.0
config.dt_init = 1.0e-9

config.max_invalid_cells = 10
config.adjust_invalid_cell_data = true
