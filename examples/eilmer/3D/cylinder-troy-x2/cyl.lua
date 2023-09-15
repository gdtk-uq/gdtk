-- cyl.lua
-- Troy and Tim's finite-length cylinder in dissociating nitrogen flow.
--
-- PJ & RJG 2016-10-16 built from the eilmer3/3D/finite-cylinder
-- and eilmer4/2D/cylinder-dlr-n90

config.dimensions = 3
D = 15.0e-3  -- diameter of cylinder, metres
L = 2.0 * D  -- (axial) length of full cylinder, will be halved later

-- Free-stream properties
T_inf = 3000.0  -- degrees K
p_inf = 2000.0  -- Pa
V_inf = 10.0e3  -- m/s
config.title = string.format("Cylinder L/D=%g in N2 at u=%g m/s.", L/D, V_inf)
print(config.title)

nsp, nmodes, gm = setGasModel('nitrogen-2sp.lua')
print("GasModel set nsp= ", nsp, " nmodes= ", nmodes)

-- Compute inflow Mach number.
Q = GasState:new{gm}
Q.p = p_inf; Q.T = T_inf; Q.massf = {N2=1.0}
gm:updateThermoFromPT(Q); gm:updateSoundSpeed(Q)
print("T=", Q.T, "density=", Q.rho, "sound speed= ", Q.a)
M_inf = V_inf / Q.a
print("M_inf=", M_inf)

inflow = FlowState:new{p=p_inf, T=T_inf, velx=V_inf, massf={N2=1.0}}
initial = FlowState:new{p=p_inf/3, T=300.0, massf={N2=1.0}}

config.reacting = true
config.reactions_file = 'e4-chem.lua'

-- Build geometry from body and flow parameters
Rc = D/2.0   -- radius of cylinder

a = Vector3:new{x=-Rc}; b = Vector3:new{y=Rc}; c = Vector3:new{x=0.0, y=0.0}

-- In order to have a grid that fits reasonably close the the shock,
-- use Billig's shock shape correlation to generate
-- a few sample points along the expected shock position.
local billig = require 'billig'
M_inf = 8.9566
print("Points on Billig's correlation.")
xys = {}
for i,y in ipairs({0.0, 0.5, 1.0, 1.5, 2.0, 2.5}) do
   x = billig.x_from_y(y*Rc, M_inf, 0.0, false, Rc)
   xys[#xys+1] = {x=x, y=y*Rc}  -- a new coordinate pair
   print("x=", x, "y=", y*Rc)
end

-- Scale the Billig distances, depending on the expected behaviour
-- relative to the gamma=1.4 ideal gas.
local b_scale = 1.1 -- for ideal (frozen-chemistry) gas
if config.reacting then
    b_scale = 0.87  -- for finite-rate chemistry
end
d = {} -- will use a list to keep the nodes for the shock boundary
for i, xy in ipairs(xys) do
   -- the outer boundary should be a little further than the shock itself
   d[#d+1] = Vector3:new{x=-b_scale*xy.x, y=b_scale*xy.y, z=0.0}
end
print("front of grid: d[1]=", d[1])

-- Extent of the cylinder in the z-direction to end face.
zshift = Vector3:new{z=L/2.0}
c2 = c + zshift
e = d[1] + zshift
f = a + zshift
g = Vector3:new{x=-Rc/2.0, y=0.0, z=L/2.0}
h = Vector3:new{x=0.0, y=Rc/2.0, z=L/2.0}
i = Vector3:new{x=0.0, y=Rc, z=L/2.0}
-- the domain is extended beyond the end of the cylinder
zshift2 = Vector3:new{z=Rc}
j = e + zshift2
k = f + zshift2

-- ...then lines, arcs, etc, that will make up the domain-end face.
xaxis = Line:new{p0=d[1], p1=a} -- first-point of shock to nose of cylinder
cylinder = Arc:new{p0=a, p1=b, centre=c}
shock = ArcLengthParameterizedPath:new{underlying_path=Spline:new{points=d}}
outlet = Line:new{p0=d[#d], p1=b}  -- top-point of shock to top of cylinder
domain_end_face = CoonsPatch:new{south=xaxis, north=outlet,
				 west=shock, east=cylinder}

-- ...lines along which we shall extrude the domain-end face
yaxis0 = Line:new{p0=d[1], p1=e}
yaxis1 = Line:new{p0=e, p1=j}

-- End-face of cylinder
xaxis = Line:new{p0=f, p1=g}
cylinder = Arc:new{p0=f, p1=i, centre=c2}
inner = Arc:new{p0=g, p1=h, centre=c2}
outlet = Line:new{p0=i, p1=h}
cyl_end_face = CoonsPatch:new{south=xaxis, north=outlet,
			      west=cylinder, east=inner}
yaxis2 = Line:new{p0=f, p1=k}

over_cylinder = SweptSurfaceVolume:new{face0123=domain_end_face,
				       edge04=yaxis0}
outside_cylinder = SweptSurfaceVolume:new{face0123=domain_end_face,
					  edge04=yaxis1}
beside_cylinder = SweptSurfaceVolume:new{face0123=cyl_end_face,
					 edge04=yaxis2}

-- Build discrete grids.
-- We choose a basic discretization and scale others from it.
nr = 20                   -- number of cells radially
nc = math.floor(1.5 * nr) -- number of cells circumferentially
na = math.floor(L/D * nc) -- number of cells along the cylinder
na1 = nc                  -- cells off the end of the cylinder
nr2 = math.floor(nr/2)    -- cells toward the cylinder axis
-- Adjust the cluster functions by trying various values.
cf0 = RobertsFunction:new{end0=true, end1=false, beta=1.5}
grid0 = StructuredGrid:new{pvolume=over_cylinder,
			   cfList={edge03=cf0, edge47=cf0},
			   niv=nr+1, njv=nc+1, nkv=na+1}
grid1 = StructuredGrid:new{pvolume=outside_cylinder,
			   cfList={edge03=cf0, edge47=cf0},
			   niv=nr+1, njv=nc+1, nkv=na1+1}
cf1 = RobertsFunction:new{end0=true, end1=false, beta=1.05}
cf2 = RobertsFunction:new{end0=true, end1=false, beta=1.6}
grid2 = StructuredGrid:new{pvolume=beside_cylinder,
			   cfList={edge01=cf1, edge32=cf2,
				   edge45=cf1, edge76=cf2},
			   niv=nr2+1, njv=nc+1, nkv=na1+1}

-- Use the grids to define some flow blocks.
-- Note that we divide up the biggest grid to make better use
-- of our multiple cpu machine.
blk0 = FBArray:new{grid=grid0, initialState=initial,
                   bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                           north=OutFlowBC_Simple:new{}},
                   nkb=math.floor(L/D)}
blk1 = FluidBlock:new{grid=grid1, initialState=initial,
		      bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
			      north=OutFlowBC_Simple:new{}}}
blk2 = FluidBlock:new{grid=grid2, initialState=initial,
		      bcList={east=OutFlowBC_Simple:new{},
			      north=OutFlowBC_Simple:new{}}}
identifyBlockConnections()

-- Set a few more config options
config.flux_calculator = "adaptive"
config.thermo_interpolator = "rhop"
config.gasdynamic_update_scheme="euler"
config.max_time = Rc/V_inf * 30
print("max_time=", config.max_time)
config.max_step = 40000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/10
-- This calculation is pretty extreme, so tolerate some bad cells.
config.adjust_invalid_cell_data = true
config.max_invalid_cells = 20
