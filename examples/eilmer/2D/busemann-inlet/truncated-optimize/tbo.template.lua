-- Author: Rowan J. Gollan and Peter J.
-- Date: 2020-09-14, 2020-09-16 multiblock
--
-- This script simulates the Busemann inlet
-- depicted in Fig. 1 in Moelder and Szpiro.
-- However, the profile of the inlet was truncated at 4 deg.
--
-- Reference:
-- Moelder and Szpiro (1966)
-- Busemann inlet for hypersonic speeds.
-- Journal of Spacecraft, 3(8):1303--1304
--
config.title = "Truncated Busemann inlet optimize trial"
config.dimensions = 2
config.axisymmetric = true
M_inf = 7.12
p_inf = 1000.0
T_inf = 350.0
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
ambient = GasState:new{gm}
ambient.p = p_inf; ambient.T = T_inf
gm:updateThermoFromPT(ambient); gm:updateSoundSpeed(ambient)
vx_inf = M_inf*ambient.a
inflow = FlowState:new{p=p_inf, T=T_inf, velx=vx_inf}
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=vx_inf}
print("vx_inf=", vx_inf)

-- Bezier points to fit truncated Busemann profile.
pts = {
   {x=3.624141, y=0.900328},
   {x=3.926667, y=$Y0}, -- 0.881750
   {x=4.452401, y=$Y1}, -- 0.826332
   {x=4.906329, y=$Y2}, -- 0.823080
   {x=5.265781, y=$Y3}, -- 0.674828
   {x=5.726706, y=$Y4}, -- 0.744892
   {x=6.159554, y=$Y5}, -- 0.508068
   {x=6.525648, y=$Y6}, -- 0.505361
   {x=6.973466, y=$Y7}, -- 0.287438
   {x=7.488979, y=$Y8}  -- 0.146429
}
contour = Bezier:new{points=pts}
wall = ArcLengthParameterizedPath:new{underlying_path=contour}
a1 = wall(0.0); a0 = {x=a1.x, y=0.0}
b1 = wall(1.0); b0 = {x=b1.x, y=0.0}
h = b1.y; isoL = 4.0*h
c1={x=b1.x+isoL, y=h}; c0={x=c1.x, y=0.0}
-- print("a1=", a1.x, a1.y)
-- print("b1=", b1.x, b1.y)
-- print("c1=", c1.x, c1.y)

quad0 = CoonsPatch:new{
   north=wall,
   east=Line:new{p0=b0, p1=b1},
   south=Line:new{p0=a0, p1=b0},
   west=Line:new{p0=a0, p1=a1}
}
quad1 = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}

grid0 = StructuredGrid:new{
   psurface=quad0, niv=201, njv=51,
   cfList={north=RobertsFunction:new{end0=false,end1=true,beta=1.2},
           south=RobertsFunction:new{end0=false,end1=true,beta=1.2}}
}
grid1 = StructuredGrid:new{
   psurface=quad1, niv=41, njv=51,
   cfList={north=RobertsFunction:new{end0=true,end1=false,beta=1.3},
           south=RobertsFunction:new{end0=true,end1=false,beta=1.3}}
}

blk0 = FBArray:new{
   grid=grid0, initialState=initial, nib=5,
   bcList={west=InFlowBC_Supersonic:new{flowState=inflow}}
}
blk1 = FluidBlock:new{
   grid=grid1, initialState=initial,
   bcList={east=OutFlowBC_Simple:new{}}
}
identifyBlockConnections()

config.max_time = 10.0e-3
config.dt_init = 1.0e-6
config.dt_plot = 1.0e-3
config.max_step = 30000
mpiDistributeBlocks{ntasks=6}
