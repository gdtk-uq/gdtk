-- cap.lua
-- Ben Stewart's Muses-C Capsule from experiments in the late 1990s.
-- PJ 2016-05-21

job_title = "X2 flow over the Muses-C aeroshell."
print(job_title)
config.dimensions = 2
config.title = job_title
config.axisymmetric = true

nsp, nmodes, gm = setGasModel('ideal-air-13-gas-model.lua')
print("GasModel set to ideal air gamma=1.3. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=4.0, T=293.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=3.3e3, T=5226.0, velx=10.3e3, vely=0.0}

-- Key points defining the geometry.
mm = 1.0e-3  -- metres per millimetre
a = Vector3:new{x=-55.0*mm, y=0.0}; b = Vector3:new{x=-55.0*mm, y=42.5*mm}
c = Vector3:new{x=-15.0*mm, y=0.0}; d = Vector3:new{x=-15.0*mm, y=15.0*mm}
e = Vector3:new{x=-12.0*mm, y=30.0*mm}; f = Vector3:new{x=0.0, y=42.5*mm}
g = Vector3:new{x=0.0, y=0.0}; h = Vector3:new{x=8.787*mm, y=21.213*mm}
i = Vector3:new{x=30.0*mm, y=0.0}; j = Vector3:new{x=17.574*mm, y=30.0*mm}
k = Vector3:new{x=30.0*mm, y=17.574*mm}; l = Vector3:new{x=30.0*mm, y=6.0*mm}
m = Vector3:new{x=42.0*mm, y=55.0*mm}; n = Vector3:new{x=45.0*mm, y=95.0*mm}
o = Vector3:new{x=45.0*mm, y=130.0*mm}; p = Vector3:new{x=0.0, y=130.0*mm}
q = Vector3:new{x=95.0*mm, y=130.0*mm}; r = Vector3:new{x=135.0*mm, y=100.0*mm}
s = Vector3:new{x=145.0*mm, y=60.0*mm}; t = Vector3:new{x=150.0*mm, y=6.0*mm}

-- Lines and arcs that run around the probe.
ab = Line:new{p0=a, p1=b}; ac = Line:new{p0=a, p1=c}
bf = Line:new{p0=b, p1=f}; cf = Bezier:new{points={c, d, e, f}};
cg = Line:new{p0=c, p1=g} 
gj = Polyline:new{segments={Arc:new{p0=g, p1=h, centre=i}, Line:new{p0=h, p1=j}}} 
fj = Line:new{p0=f, p1=j}
fp = Line:new{p0=f, p1=p}; po = Line:new{p0=p, p1=o}
jo = Bezier:new{points={j, m, n, o}}
jk = Line:new{p0=j, p1=k}
oos = Bezier:new{points={o, q, r, s}}
ks = Line:new{p0=k, p1=s}; kl = Line:new{p0=k, p1=l}
lt = Line:new{p0=l, p1=t}; st = Line:new{p0=s, p1=t}

patch0 = makePatch{north=bf, east=cf, south=ac, west=ab, gridType="ao"}
patch1 = makePatch{north=fj, east=gj, south=cg, west=cf, gridType="ao"}
patch2 = makePatch{north=po, east=jo, south=fj, west=fp, gridType="ao"}
patch3 = makePatch{north=oos, east=ks, south=jk, west=jo, gridType="tfi"}
patch4 = makePatch{north=st, east=lt, south=kl, west=ks, gridType="tfi"}

-- Mesh the patches, with particular discretisation.
-- Some fiddling with the clustering functions is required
-- to get the cells to distribute nicely in front of the aeroshell.
-- Don't need to be so fussy downstream and away from the body.
factor = 2
nx0 = 40*factor; ny0 = 60*factor
cfxn0 = RobertsFunction:new{end0=false,end1=true,beta=1.4}
cfxs0 = RobertsFunction:new{end0=false,end1=true,beta=1.2}
grid0 = StructuredGrid:new{psurface=patch0, niv=nx0+1, njv=ny0+1,
			   cfList={north=cfxn0, south=cfxs0}}
nx1 = 30*factor
grid1 = StructuredGrid:new{psurface=patch1, niv=nx1+1, njv=ny0+1}
ny2 = 60*factor
cfy2 = RobertsFunction:new{end0=true,end1=false,beta=1.2}
grid2 = StructuredGrid:new{psurface=patch2, niv=nx1+1, njv=ny2+1,
			   cfList={east=cfy2, west=cfy2}}
grid3 = StructuredGrid:new{psurface=patch3, niv=nx1+1, njv=ny2+1,
			   cfList={east=cfy2, west=cfy2}}
grid4 = StructuredGrid:new{psurface=patch4, niv=nx1+1, njv=ny2+1,
			   cfList={east=cfy2, west=cfy2}}

-- Define the flow-solution blocks.
blk0 = FluidBlock:new{grid=grid0, fillCondition=initial, label="BLK0"}
blk1 = FluidBlock:new{grid=grid1, fillCondition=initial, label="BLK1"}
blk2 = FluidBlock:new{grid=grid2, fillCondition=initial, label="BLK2"}
blk3 = FluidBlock:new{grid=grid3, fillCondition=initial, label="BLK3"}
blk4 = FluidBlock:new{grid=grid4, fillCondition=initial, label="BLK4"}
-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_Supersonic:new{flowCondition=inflow, label="inflow-boundary"}
blk2.bcList[north] = OutFlowBC_Simple:new{label="outflow-boundary"}
blk3.bcList[north] = OutFlowBC_Simple:new{label="outflow-boundary"}
blk4.bcList[north] = OutFlowBC_Simple:new{label="outflow-boundary"}

-- Do a little more setting of global data.
config.max_time = 80.0e-6  -- seconds
config.max_step = 50000
config.dt_init = 0.5e-8
config.cfl_value = 0.5
config.dt_plot = 5.0e-6
