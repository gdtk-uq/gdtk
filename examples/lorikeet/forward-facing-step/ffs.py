# ffs.py
# Supersonic flow over a forward-facing step.
# PJ 2022-10-08 adapted from the Eilmer example.
#    2022-12-11 back to being a pure 2D example for Lorikeet.
#
config.title = "Forward-facing step with supersonic flow."
print(config.title)
#
init_gas_model("ideal-air-gas-model.lua")
initial = FlowState(p=101.325e3, T=300.0)
sspeed = initial.gas.a
print("sound speed=", sspeed)
inflow = FlowState(p=101.325e3, T=300.0, velx=3.0*sspeed)
#
# Plan-view of the flow domain. k,z is out of the page.
#
#  ^j,y
#  |
#  a2-----b2---------------c2
#  |      |                |
#  |   1  |        2       |
#  |      |                |
#  a1-----b1---------------c1
#  |   0  |
#  a0-----b0   ---->i,x
#
a0 = Vector3(0.0, 0.0); a1 = Vector3(0.0, 0.2); a2 = Vector3(0.0, 1.0)
b0 = Vector3(0.6, 0.0); b1 = Vector3(0.6, 0.2); b2 = Vector3(0.6, 1.0)
c0 = Vector3(3.0, 0.0); c1 = Vector3(3.0, 0.2); c2 = Vector3(3.0, 1.0)
#
box0 = CoonsPatch(p00=a0, p10=b0, p11=b1, p01=a1)
box1 = CoonsPatch(p00=a1, p10=b1, p11=b2, p01=a2)
box2 = CoonsPatch(p00=b1, p10=c1, p11=c2, p01=b2)
#
# Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = int(0.6/dx); nbc = int(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = int(0.2/dx); n12 = int(0.8/dx)
print("n01=", n01, "n12=", n12)
#
grd0 = StructuredGrid(psurf=box0, niv=nab+1, njv=n01+1)
grd1 = StructuredGrid(psurf=box1, niv=nab+1, njv=n12+1)
grd2 = StructuredGrid(psurf=box2, niv=nbc+1, njv=n12+1)
#
blk0 = FluidBlock(i=0, j=0, grid=grd0, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow),'jplus':ExchangeBC()})
blk1 = makeFBArray(i0=0, j0=1, nj=4, grid=grd1, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow), 'jminus':ExchangeBC(),'iplus':ExchangeBC()})
blk2 = makeFBArray(i0=1, j0=1, ni=4, nj=4, grid=grd2, initialState=inflow,
                   bcs={'iminus':ExchangeBC(),'iplus':OutflowBC()})
#
# Although the sumulation code can cope with not having a full array of blocks,
# the StructuredGrid in Paraview cannot tolerate missing pieces.
# To deal with that, we will make a dummy volume and carry it along.
# Once we get the post-processor to build and unstructured-grid for Paraview,
# we can drop it.
box3 = CoonsPatch(p00=b0, p10=c0, p11=c1, p01=b1)
grd3 = StructuredGrid(psurf=box3, niv=nbc+1, njv=n01+1)
blk3 = makeFBArray(i0=1, j0=0, ni=4, grid=grd3, initialState=initial, active=False)
#
config.max_time = 5.0e-3
config.max_step = 6000
