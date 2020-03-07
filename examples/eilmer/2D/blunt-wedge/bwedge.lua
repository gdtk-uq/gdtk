-- bwedge.lua
-- To model the flow shown in Figure 10 of Mike Macrossan's paper
-- "Hypervelocity flow of dissociating nitrogen downstream of a blunt nose"
-- Journal of Fluid Mechanics (1990) Vol. 217, pages 167--202.
-- PJ 2020-03-07
--
config.title = "Dissociating nitrogen flow over a blunt wedge."
print(config.title)

nsp, nmodes, gm = setGasModel('nitrogen-2sp.lua')
print("GasModel set nsp=", nsp, "nmodes=", nmodes)
config.reactions_file = 'n2-n-reactions.lua'
config.reacting = true

gs = GasState:new{gm} -- Note the single table argument!
gs.rho = 4.41e-2 -- kg/m^3
gs.T = 4415 -- degrees K
molef = {N=0.172, N2=1-0.172} -- dissociation fraction 0.094
gs.massf = gm:molef2massf(molef)
gm:updateThermoFromRHOT(gs)
gm:updateSoundSpeed(gs)
print("Free stream p=", gs.p, "T=", gs.T, "massf-N2=", gs.massf["N2"])
inflow = FlowState:new{p=gs.p, T=gs.T, velx=6360.0, massf=gs.massf}
initial = FlowState:new{p=5.0, T=300.0, massf={N2=1.0, N=0.0}}
Minf = inflow.velx / inflow.a
print("Minf=", Minf)

-- First, we will define a patch around the nose of the plate.
Rn = 0.005 -- radius of cylindrical nose, m
xEnd = 10*Rn -- downstream extent of wedge
alpha = math.rad(15) -- angle of wedge wrt free stream
delta = 0.002 -- offset for inflow boundary, to accommodate detached shock
-- Second, specify surface of cylinder and wedge.
a = {x=0, y=0} -- Centre of curvature for nose
b = {x=-Rn, y=0}
c = {x=-Rn*math.sin(alpha), y=Rn*math.cos(alpha)}
bc = Arc:new{p0=b, p1=c, centre=a}
--  Down-stream end of wedge
d = {x=xEnd, y=c.y+(xEnd-c.x)*math.tan(alpha)}
cd = Line:new{p0=c, p1=d}
-- Outer-edge of flow domain has to contain the shock layer.
-- Allow sufficient for shock stand-off at the stagnation line.
R2 = Rn + delta
e = {x=-R2, y=0}
-- The shock angle is eventually fairly low but the blunt nose
-- displaces the shock a long way out, so we allow some more space.
-- We need to set the boundary high enough to avoid the shock.
R3 = 2*Rn
f = {x=-R3*math.sin(alpha), y=R3*math.cos(alpha)}
-- Now, put in intermediate control points so that we can use
-- cubic Bezier curve for the inflow boundary around the nose
-- and a straight line downstream of point f.
e1 = {x=e.x, y=Rn}
alpha2 = math.rad(35) -- guessed angle for straight boundary
e2 = {x=f.x-delta*math.cos(alpha2), y=f.y-delta*math.sin(alpha2)}
ef = Bezier:new{points={e, e1, e2, f}}
g = {x=xEnd, y=f.y+(xEnd-f.x)*math.tan(alpha2)}
fg = Line:new{p0=f, p1=g}
-- Define straight-line segments between surface and outer boundary.
eb = Line:new{p0=e, p1=b}
fc = Line:new{p0=f, p1=c}; cf = ReversedPath:new{underlying_path=fc}
dg = Line:new{p0=d, p1=g}
--
-- Define the blocks using the path segments.
-- Note that the east face of region0 wraps around the nose and
-- that the north face of region0 is adjacent to the west face
-- of region1.
region0 = CoonsPatch:new{north=fc, east=bc, south=eb, west=ef}
region1 = CoonsPatch:new{north=fg, east=dg, south=cd, west=cf}
cfun0 = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grid0 = StructuredGrid:new{psurface=region0, niv=41, njv=41,
                           cfList={north=cfun0, south=cfun0}}
cfun1 = RobertsFunction:new{end0=true, end1=false, beta=1.2}
cfun2 = RobertsFunction:new{end0=true, end1=false, beta=1.2}
grid1 = StructuredGrid:new{psurface=region1, niv=101, njv=41,
                           cfList={east=cfun1, west=cfun1,
                                   north=cfun2, south=cfun2}}
fba0 = FBArray:new{grid=grid0, initialState=initial, label="nose",
                   bcList={west=InFlowBC_Supersonic:new{flowState=inflow}}, 
                   nib=1, njb=2}
fba1 = FBArray:new{grid=grid1, initialState=initial, label="wedge",
                   bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
                           east=OutFlowBC_Simple:new{}}, 
                   nib=5, njb=1}
identifyBlockConnections()

-- Set a few more config options.
-- To get a reasonable start, we needed to set dt_init.
config.max_time = 100.0e-6
config.max_step = 40000
config.dt_plot = 20.0e-6
config.dt_init = 5.0e-9
