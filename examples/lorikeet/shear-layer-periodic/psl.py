# psl.py
# Periodic shear layer.
# PJ 2022-10-09 Adapted from Eilmer4 example.
#    2022-12-15 Lorikeet variant.
#
config.title = "Periodic shear layer."
print(config.title)
init_gas_model("ideal-air-gas-model.lua")
#
H = 0.010 # y-direction layer thickness in metres
L = 0.100 # x-direction wavelength in metres
ymin = -20.0*H; ymax = 20.0*H
xmin = -L; xmax = L
#
box0 = CoonsPatch(p00=Vector3(xmin,ymin), p10=Vector3(xmax,ymin),
                  p11=Vector3(xmax,ymax), p01=Vector3(xmin,ymax))
factor = 1
grd0 = StructuredGrid(psurf=box0, niv=int(60*factor)+1, njv=int(120*factor)+1)
#
def initial_flow(x, y):
    """
    User-defined function for the initial flow state works in physical space.
    """
    p = 100.0e3 # Pa
    T = 300.0 # K
    velx0 = 200.0 # m/s subsonic
    # The lower half of the domain is flowing left and the upper half, right.
    velx = -velx0
    if y > H:
        velx = velx0
    elif y > -H:
        velx = y/H * velx0
    # Add perturbation that is periodic west to east
    # but gets smaller toward the north and south boundaries.
    vely = 10.0 * math.exp(-abs(y)/H) * \
        (math.cos(x/L*math.pi) + math.sin(2*x/L*math.pi))
    # We use the FlowState object to conveniently set all of
    # the relevant properties.
    return FlowState(p=p, velx=velx, vely=vely, T=T)
#
blk0 = FluidBlock(grid=grd0, initialState=initial_flow,
                  bcs={'iminus':ExchangeBC(),'iplus':ExchangeBC()})
#
config.max_time = 10.0e-3  # seconds
config.max_step = 150000
config.plot_dt = 0.1e-3
