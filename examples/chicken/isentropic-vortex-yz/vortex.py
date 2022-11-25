# vortex.py
# Isentropic Vortex.
# CM 2022-10-18 Adapted from Eilmer4 example.
#
config.title = "Isentropic Vortex."
print(config.title)
#
H = 10.0 # y-direction in metres
L = 10.0 # x-direction in metres
ymin = -0.2*L; ymax = 1.2*L
xmin = -0.2; xmax = 0.0
zmin = -0.2*H; zmax = 1.2*H
#
vol0 = TFIVolume(p000=Vector3(xmin,ymin,zmin), p100=Vector3(xmax,ymin,zmin),
                 p110=Vector3(xmax,ymax,zmin), p010=Vector3(xmin,ymax,zmin),
                 p001=Vector3(xmin,ymin,zmax), p101=Vector3(xmax,ymin,zmax),
                 p111=Vector3(xmax,ymax,zmax), p011=Vector3(xmin,ymax,zmax))
factor = 1
grd0 = StructuredGrid(pvolume=vol0, niv=3, njv=int(120*factor)+1, nkv=int(120*factor)+1)
#
def initial_flow(x, y, z):
    """
    User-defined function for the initial flow state works in physical space.
    """
    # Ideal air gas model
    Rgas = 287.1; Gamma = 1.4
    # Undisturbed gas state
    p0 = 100.0e3 # Pa
    T0 = 300.0 # K
    rho0 = p0/(Rgas*T0)
    a0 = math.sqrt(Gamma*Rgas*T0)
    # Vortex centre and strength
    yc = L/2.0; zc = H/2.0
    eps = 5.0*a0
    ybar = y - yc; zbar = z - zc
    r2 = ybar**2 + zbar**2
    r = math.sqrt(r2)
    # Vortex flow
    tmp0 = (eps/(2*math.pi))**2 * math.exp(1.0-r2)
    v2 = tmp0*r2
    v = math.sqrt(v2)
    dT = (Gamma-1.0)/(Rgas*Gamma)*0.5*tmp0
    T = T0 - dT
    tmp1 = math.pow(rho0, Gamma)/p0 * Rgas * T
    rho = math.pow(tmp1, 1.0/(Gamma-1.0))
    pressure = rho*Rgas*T
    if r < 1.0e-9:
        vely = 0.0; velz = 0.0
    else:
        vely = -zbar/r * v; velz = ybar/r * v
    # Add background velocity
    vely += a0; velz += a0
    return FlowState(p=pressure, velx=0.0, vely=vely, velz=velz, T=T)
#
blk0 = FluidBlock(grid=grd0, initialState=initial_flow,
                  bcs={'jminus':ExchangeBC(),'jplus':ExchangeBC(),
                       'kminus':ExchangeBC(),'kplus':ExchangeBC()})
#
config.max_time = 50.0e-3  # seconds
config.max_step = 150000
config.flux_calc = "sbp_asf"
add_dt_plot(0.0, 1.0e-3)
