# ffs.py
# Supersonic flow over a forward-facing step.
# PJ 2022-10-08 adapted from the Eilmer example.
#
config.title = "Forward-facing step with supersonic flow."
print(config.title)
#
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
a0 = Vector3(0.0, 0.0, 0.0); a1 = Vector3(0.0, 0.2, 0.0); a2 = Vector3(0.0, 1.0, 0.0)
b0 = Vector3(0.6, 0.0, 0.0); b1 = Vector3(0.6, 0.2, 0.0); b2 = Vector3(0.6, 1.0, 0.0)
c0 = Vector3(3.0, 0.0, 0.0); c1 = Vector3(3.0, 0.2, 0.0); c2 = Vector3(3.0, 1.0, 0.0)
H = Vector3(0.0, 0.0, 0.1)
vol0 = TFIVolume(p000=a0, p100=b0, p110=b1, p010=a1,
                 p001=a0+H, p101=b0+H, p111=b1+H, p011=a1+H)
vol1 = TFIVolume(p000=a1, p100=b1, p110=b2, p010=a2,
                 p001=a1+H, p101=b1+H, p111=b2+H, p011=a2+H)
vol2 = TFIVolume(p000=b1, p100=c1, p110=c2, p010=b2,
                 p001=b1+H, p101=c1+H, p111=c2+H, p011=b2+H)
#
# Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = int(0.6/dx); nbc = int(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = int(0.2/dx); n12 = int(0.8/dx)
print("n01=", n01, "n12=", n12)
#
grd0 = StructuredGrid(pvolume=vol0, niv=nab+1, njv=n01+1, nkv=3)
grd1 = StructuredGrid(pvolume=vol1, niv=nab+1, njv=n12+1, nkv=3)
grd2 = StructuredGrid(pvolume=vol2, niv=nbc+1, njv=n12+1, nkv=3)
#
blk0 = FluidBlock(i=0, j=0, grid=grd0, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow),'jplus':ExchangeBC()})
blk1 = FluidBlock(i=0, j=1, grid=grd1, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow), 'jminus':ExchangeBC(),'iplus':ExchangeBC()})
blk2 = FluidBlock(i=1, j=1, grid=grd2, initialState=inflow,
                  bcs={'iminus':ExchangeBC(),'iplus':OutflowBC()})
#
# Although the sumulation code can cope with not having a full array of blocks,
# the StructuredGrid in Paraview cannot tolerate missing pieces.
# To deal with that, we will make a dummy volume and carry it along.
# Once we get the post-processor to build and unstructured-grid for Paraview,
# we can drop it.
vol3 = TFIVolume(p000=b0, p100=c0, p110=c1, p010=b1,
                 p001=b0+H, p101=c0+H, p111=c1+H, p011=b1+H)
grd3 = StructuredGrid(pvolume=vol3, niv=nbc+1, njv=n01+1, nkv=3)
blk3 = FluidBlock(i=1, j=0, grid=grd3, initialState=initial, active=False)
#
config.max_time = 5.0e-3
config.max_step = 6000
add_cfl_value(0.0, 0.5)
add_dt_plot(0.0, 1.0e-3)

