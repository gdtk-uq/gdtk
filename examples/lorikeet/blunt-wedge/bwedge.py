# bwedge.lua
# To model the flow shown in Figure 10 of Mike Macrossan's paper
# "Hypervelocity flow of dissociating nitrogen downstream of a blunt nose"
# Journal of Fluid Mechanics (1990) Vol. 217, pages 167--202.
# PJ 2020-03-07
#    2023-01-01 Adapted from the Eilmer4 example.
#
config.title = "Dissociating nitrogen flow over a blunt wedge."
print(config.title)
#
init_gas_model('nitrogen-2sp.lua')
gs = GasState(config.gmodel)
gs.rho = 4.41e-2 # kg/m^3
gs.T = 4415 # degrees K
gs.molef = {'N':0.172, 'N2':1-0.172} # dissociation fraction 0.094
gs.update_thermo_from_rhoT()
gs.update_sound_speed()
print("Free stream p=", gs.p, "T=", gs.T, "massf-N2=", gs.massf_as_dict["N2"])
inflow = FlowState(p=gs.p, T=gs.T, velx=6360.0, massf=gs.massf)
initial = FlowState(p=5.0, T=300.0, velx=0.0, massf={'N2':1.0, 'N':0.0})
Minf = inflow.vel.x / inflow.gas.a
print("Minf=", Minf)
config.reacting = True
config.reaction_file_1 = 'n2-n-reactions.lua'
#
# Geometry
from gdtk.geom.path import Arc, Line, Bezier
from gdtk.geom.cluster import RobertsFunction
from gdtk.geom.surface import CoonsPatch
#
# First, we will define a patch around the nose of the plate.
Rn = 0.005 # radius of cylindrical nose, m
xEnd = 10*Rn # downstream extent of wedge
alpha = math.radians(15) # angle of wedge wrt free stream
delta = 0.002 # offset for inflow boundary, to accommodate detached shock
# Second, specify surface of cylinder and wedge.
c = Vector3(x=0, y=0) # Centre of curvature for nose
a = Vector3(x=-Rn, y=0)
b = Vector3(x=-Rn*math.sin(alpha), y=Rn*math.cos(alpha))
ab = Arc(a=a, b=b, c=c)
#  Down-stream end of wedge
d = Vector3(x=xEnd, y=b.y+(xEnd-b.x)*math.tan(alpha))
bd = Line(p0=b, p1=d)
# Outer-edge of flow domain has to contain the shock layer.
# Allow sufficient for shock stand-off at the stagnation line.
R2 = Rn + delta
e = Vector3(x=-R2, y=0)
# The shock angle is eventually fairly low but the blunt nose
# displaces the shock a long way out, so we allow some more space.
# We need to set the boundary high enough to avoid the shock.
R3 = 2*Rn
f = Vector3(x=-R3*math.sin(alpha), y=R3*math.cos(alpha))
# Now, put in intermediate control points so that we can use
# cubic Bezier curve for the inflow boundary around the nose
# and a straight line downstream of point f.
e1 = Vector3(x=e.x, y=Rn)
alpha2 = math.radians(35) # guessed angle for straight boundary
e2 = Vector3(x=f.x-delta*math.cos(alpha2), y=f.y-delta*math.sin(alpha2))
ef = Bezier(B=[e, e1, e2, f])
g = Vector3(x=xEnd, y=f.y+(xEnd-f.x)*math.tan(alpha2))
fg = Line(p0=f, p1=g)
# Define straight-line segments between surface and outer boundary.
ae = Line(p0=a, p1=e)
bf = Line(p0=b, p1=f)
dg = Line(p0=d, p1=g)
#
# Define the blocks using the path segments.
# Note that the south face of region0 wraps around the nose and that
# the east face of region0 is adjacent to the west face of region1.
region0 = CoonsPatch(north=ef, east=bf, south=ab, west=ae)
region1 = CoonsPatch(north=fg, east=dg, south=bd, west=bf)
cfun0 = RobertsFunction(end0=True, end1=False, beta=1.2)
grid0 = StructuredGrid(psurf=region0, niv=41, njv=41,
                       cf_list=[cfun0, None, None, None])
cfuni = RobertsFunction(end0=True, end1=False, beta=1.2)
grid1 = StructuredGrid(psurf=region1, niv=101, njv=41,
                       cf_list=[cfuni, None, cfuni, None])
#
blk0 = FluidBlock(i=0, grid=grid0, initialState=initial,
                  bcs={'north':InflowBC(inflow),'east':ExchangeBC()})
blk1 = FluidBlock(i=1, grid=grid1, initialState=initial,
                  bcs={'west':ExchangeBC(),'north':InflowBC(inflow), 'east':OutflowBC()})
#
config.flux_calc = FluxCalc.ausmdv_plus_hanel
config.max_time = 100.0e-6
config.max_step = 40000
config.plot_dt = 20.0e-6
# To get a reasonable start, we needed to set dt_init.
config.dt_init = 1.0e-9
