# sod.py
# Simulate the Sod shock tube with chicken.
# PJ 2022-09-19
# We'll model the shock tube as a pair of simple boxes.
H = 1.0 # y-direction in metres
L = 1.0 # x-direction in metres
ymin = -0.2; ymax = 0.0
xmin = 0.0; xmax = L
zmin = 0.0 ; zmax = H
vol0 = TFIVolume(p000=Vector3(xmin,ymin,zmin), p100=Vector3(xmax,ymin,zmin),
                 p110=Vector3(xmax,ymax,zmin), p010=Vector3(xmin,ymax,zmin),
                 p001=Vector3(xmin,ymin,zmax), p101=Vector3(xmax,ymin,zmax),
                 p111=Vector3(xmax,ymax,zmax), p011=Vector3(xmin,ymax,zmax))
factor = 1
grd0 = StructuredGrid(pvolume=vol0, niv=int(120*factor)+1, njv=3, nkv=int(120*factor)+1)
#
def initial_flow(x, y, z):
    high = FlowState(p=100.0e3, T=348.4)
    low = FlowState(p=10.0e3, T=278.8)
    r = math.sqrt(math.pow(x-0.5,2)+math.pow(z-0.5,2))
    if r < 0.25:
        return high
    else:
        return low       
#
blk0 = FluidBlock(grid=grd0, initialState=initial_flow)
config.max_time = 0.6e-3
config.max_step = 1000
