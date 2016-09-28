-- ram2.lua
-- Ramjet design from Birte Haker, with ReactionZones
-- Eilmer versions PJ, 25-Mar-2007, 02-Jun-2007, 22-Apr-2011
-- Eilmer4 version PJ, 26-Sep-2016
--

function radians(degrees)
   return degrees/180.0*math.pi
end
function polar_vector(degrees)
   return Vector3:new{x=math.cos(radians(degrees)),y=math.sin(radians(degrees))}
end

config.title = "Ramjet in Mach 1.8 Flight -- Reacting Gas."
print(config.title)

-- Gas model setup taken from Bittker combsution in a duct example.
nsp, nmodes, gm = setGasModel("h2-o2-n2-9sp.lua")
config.reacting = true
config.reactions_file = "h2-o2-n2-9sp-18r.lua"
molefInit = {O2=0.1480, N2=0.5562, H2=0.2958}
massfInit = gm:molef2massf(molefInit)
print("GasModel set to H2-O2-N2 mix. nsp= ", nsp, " nmodes= ", nmodes)

-- Define flow conditions
p_inf = 90.0e3  -- Pa
T_inf = 282     -- degrees K
a_inf = math.sqrt(1.4 * 287.0 * T_inf) -- assumes ideal air
M_inf = 1.8
Vx_inf = M_inf * a_inf
print("flight speed=", Vx_inf, "m/s")
initial = FlowState:new{p=p_inf/3, T=T_inf, velx=0.0, massf=massfInit}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=Vx_inf, massf=massfInit}

-- Set up the geometry in the (x,y)-plane be first defining key parameters,
-- then defining points and paths using those parameters.
mm = 1.0e-3         -- metres per mm
L0 = 70.0 * mm      -- distance upstream of inlet
Y0 = 300.0 * mm     -- outer edge of flow domain
Y1 = 55.0/2 * mm    -- inlet
Y2 = 63.0 * mm      -- maximum width of body
L1 = 120.0 * mm     -- forebody length
L2 = 300.0 * mm     -- combustor length
L3 = 30.0 * mm      -- nozzle length
L4 = 100.0 * mm     -- trailing body
-- Center-body
R2 = 24.0 * mm      -- bluff nose
R3 = R2 + 20.0 * mm -- block thickness around nose
Y3 = 21.0 * mm      -- centrebody diameter

-- Derived parameters.
-- R1 = raduis of curvature of forebody
R1 = (L1 * L1 + (Y2 - Y1)^2) / (2.0 * (Y2 - Y1))
theta1 = math.asin(L1 / R1)
print("R1=", R1, "theta1=", theta1)

-- Key points on (x,y)-plane
-- Outside of cowl
a = Vector3:new{x=0.0, y=Y1}
b = Vector3:new{x=L1, y=Y2}
c = Vector3:new{x=L1, y=Y2-R1} -- centre of curvature
z0 = Vector3:new{x=-L0, y=0.0}; z1 = Vector3:new{x=-L0, y=2*Y1/3.0}
z2 = Vector3:new{x=-L0, y=Y1}; z3 = Vector3:new{x=-L0, y=Y0}
z4 = Vector3:new{x=0.0, y=0.0}; z5 = Vector3:new{x=0.0, y=2*Y1/3.0}
z6 = Vector3:new{x=0.0, y=Y0}
z7 = Vector3:new{x=L1, y=Y0}; z8 = Vector3:new{x=L1+L2+L3, y=Y0}
z9 = z8 + Vector3:new{x=L4, y=0.0}
b2 = b + Vector3:new{x=L2+L3, y=0.0}
z10 = b2 + Vector3:new{x=L4, y=0.0}
-- Inside of cowl.
g = b - Vector3:new{x=0.0, y=5.0*mm}
h = 0.33*a + 0.67*g  -- part way along the line
f = g + Vector3:new{x=L2-5.0*mm, y=0.0}
e = f + Vector3:new{x=0.0, y=-5.0*mm}
-- d is defined after q

-- Centrebody
i = Vector3:new{x=L1-5.0*mm, y=0.0}
j = i - Vector3:new{x=R2, y=0.0}
k = i + R2 * polar_vector(135.0)
l = i + Vector3:new{x=0.0, y=R2}
m = Vector3:new{x=L1, y=R2}
n = Vector3:new{x=L1, y=Y3}
o = Vector3:new{x=L1+L2, y=Y3}
p = o + Vector3:new{x=0.0, y=5.0*mm}
q = p + 20.0 * mm * polar_vector(30.0)
r = q + 5.0 * mm * polar_vector(10.0)
s = q + 20.0 * mm * polar_vector(5.0)
t = Vector3:new{x=z10.x, y=s.y}

-- Finish throat of nozzle
d = Vector3:new{x=q.x, y=e.y}

-- Intermediate points through combustor
z11 = i - Vector3:new{x=R3, y=0.0}
z12 = i + R3 * polar_vector(135.0)
z13 = i + Vector3:new{x=0.0, y=R3}
z14 = Vector3:new{x=L1, y=R3}

quad = {}
-- Region upstream of inlet.
quad[0] = CoonsPatch:new{p11=z5, p10=z4, p00=z0, p01=z1}
quad[1] = CoonsPatch:new{p11=a, p10=z5, p00=z1, p01=z2}
quad[2] = CoonsPatch:new{p11=z6, p10=a, p00=z2, p01=z3}

-- Region around the outer cowl.
quad[3] = CoonsPatch:new{north=Line:new{p0=z6, p1=z7},
			 east=Line:new{p0=b, p1=z7},
			 south=Arc:new{p0=a, p1=b, centre=c},
			 west=Line:new{p0=a, p1=z6}}
quad[4] = CoonsPatch:new{p11=z8, p10=b2, p00=b, p01=z7}

-- Region inside entrance
quad[5] = CoonsPatch:new{north=Line:new{p0=z5, p1=z12},
			 east=Arc:new{p0=z11, p1=z12, centre=i},
			 south=Line:new{p0=z4, p1=z11},
			 west=Line:new{p0=z4, p1=z5}}
quad[6] = CoonsPatch:new{p11=h, p10=z12, p00=z5, p01=a}
quad[7] = CoonsPatch:new{north=Arc:new{p0=z11, p1=z12, centre=i},
			 east=Line:new{p0=k, p1=z12},
			 south=Arc:new{p0=j, p1=k, centre=i},
			 west=Line:new{p0=j, p1=z11}}
south8 = Polyline:new{segments={Arc:new{p0=k, p1=l, centre=i},
				Line:new{p0=l, p1=m}}}
south9north8 = Polyline:new{segments={Arc:new{p0=z12, p1=z13, centre=i},
				      Line:new{p0=z13, p1=z14}}}
quad[8] = CoonsPatch:new{north=south9north8,
			 east=Line:new{p0=m, p1=z14},
			 south=south8,
			 west=Line:new{p0=k, p1=z12}}
quad[9] = CoonsPatch:new{north=Line:new{p0=h, p1=g},
			 east=Line:new{p0=z14, p1=g},
			 south=south9north8,
			 west=Line:new{p0=z12, p1=h}}

-- Long part of combustor
quad[10] = CoonsPatch:new{p11=p, p10=o, p00=n, p01=m}
quad[11] = CoonsPatch:new{p11=e, p10=p, p00=m, p01=z14}
quad[12] = CoonsPatch:new{p11=f, p10=e, p00=z14, p01=g}

-- Converging-diverging nozzle
north13 = Polyline:new{segments={Line:new{p0=e, p1=d},
				 Line:new{p0=d, p1=b2}}}
south13 = Polyline:new{segments={Line:new{p0=p, p1=q},
				 Line:new{p0=q, p1=r},
				 Line:new{p0=r, p1=s}}}
quad[13] = CoonsPatch:new{north=north13,
			  east=Line:new{p0=s, p1=b2},
			  south=south13,
			  west=Line:new{p0=p, p1=e}}

-- Blocks after nozzle.
quad[14] = CoonsPatch:new{p11=z10, p10=t, p00=s, p01=b2}
quad[15] = CoonsPatch:new{p11=z9, p10=z10, p00=b2, p01=z8}

-- Discretise the patches.
nx0 = 30; ny0 = 20; ny1 = 10; ny2 = 70
nx1 = 50; nx2 = math.floor(nx1*0.6); nx3 = 4*nx1; ny3 = 20; ny4 = 4
nx4 = 40

grid = {}
grid[0] = StructuredGrid:new{psurface=quad[0], niv=nx0+1, njv=ny0+1}
grid[1] = StructuredGrid:new{psurface=quad[1], niv=nx0+1, njv=ny1+1}
grid[2] = StructuredGrid:new{psurface=quad[2], niv=nx0+1, njv=ny2+1}
grid[3] = StructuredGrid:new{psurface=quad[3], niv=nx2+nx2+1, njv=ny2+1}
grid[4] = StructuredGrid:new{psurface=quad[4], niv=nx3+1, njv=ny2+1}
grid[5] = StructuredGrid:new{psurface=quad[5], niv=nx1+1, njv=ny0+1}
grid[6] = StructuredGrid:new{psurface=quad[6], niv=nx1+1, njv=ny1+1}
grid[7] = StructuredGrid:new{psurface=quad[7], niv=ny0+1, njv=ny3+1}
grid[8] = StructuredGrid:new{psurface=quad[8], niv=nx2+1, njv=ny3+1}
grid[9] = StructuredGrid:new{psurface=quad[9], niv=nx2+1, njv=ny1+1}
grid[10] = StructuredGrid:new{psurface=quad[10], niv=nx3+1, njv=ny4+1}
grid[11] = StructuredGrid:new{psurface=quad[11], niv=nx3+1, njv=ny3+1}
grid[12] = StructuredGrid:new{psurface=quad[12], niv=nx3+1, njv=ny1+1}
grid[13] = StructuredGrid:new{psurface=quad[13], niv=nx4+1, njv=ny3+1}
grid[14] = StructuredGrid:new{psurface=quad[14], niv=nx3+1, njv=ny3+1}
grid[15] = StructuredGrid:new{psurface=quad[15], niv=nx3+1, njv=ny2+1}

blk = {}
for ib = 0, 15 do
   blk[ib] = SBlock:new{grid=grid[ib], fillCondition=initial}
end

-- Boundary conditions
identifyBlockConnections()
for ib = 0, 2 do
   blk[ib].bcList[west] = InFlowBC_Supersonic:new{flowCondition=inflow}
end
blk[14].bcList[east] = OutFlowBC_Simple:new{}
blk[15].bcList[east] = OutFlowBC_Simple:new{}

-- Do a little more setting of global data.
config.axisymmetric = true
config.flux_calculator = "adaptive"
config.max_time = 5.0e-3  -- seconds
print("maximum simulation time=", config.max_time, "seconds")
config.max_step = 80000
config.dt_init = 1.0e-8
config.dt_plot = config.max_time / 20.0  -- so that we get a few images
config.dt_history = 10.0e-5

-- We'll let reaction occur in the downstream end of the combustor
-- and into the nozzle.
ReactionZone:new{p0=Vector3:new{x=L1+L2/4.0, y=Y3},
		 p1=Vector3:new{x=L1+L2+L3, y=Y2}}
IgnitionZone:new{p0=Vector3:new{x=L1+L2/4.0, y=Y3},
		 p1=Vector3:new{x=L1+L2/3.0, y=Y2}, T=1100.0}
-- Also, we want the flow through the combustor to establish before
-- allowing reactions to become active.
config.reaction_time_delay = 2.5e-3

setHistoryPoint{x=10.0*mm, y=2.0*mm}  -- just inside the inlet

dofile("sketch-domain.lua")
