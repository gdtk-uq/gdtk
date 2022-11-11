# kh.py
# A simple 2D Kelvin-Helmholtz instability.
# PJ 2022-10-09 Adapted from Lachlan Whyborn's Eilmer4 example.
#
import math
config.title = "Kelvin-Helmholtz instability."
print(config.title)
#
H = 0.010 # z-direction layer thickness in metres
L = 1.000 # x,y-direction domain size in metres
#
vol0 = TFIVolume(p000=Vector3(0,0,0), p100=Vector3(L,0,0),
                 p110=Vector3(L,L,0), p010=Vector3(0,L,0),
                 p001=Vector3(0,0,H), p101=Vector3(L,0,H),
                 p111=Vector3(L,L,H), p011=Vector3(0,L,H))
N = 256
grd0 = StructuredGrid(pvolume=vol0, niv=N+1, njv=N+1, nkv=3)
#
# Flow conditions
U1 = 50.0; U2 = -50.0; U_m = (U1-U2)/2
rho1 = 1.0; rho2 = 2.0; rho_m = (rho1-rho2)/2
p_inf = 1.0e5
Rgas = 287.1; g = 1.4
Ls = 0.025

def initial_flow(x, y, z):
    """
    User-defined function for the initial flow state works in physical space.
    """
    if y < L/4:
        rho = rho1 - rho_m * math.exp((y-L/4)/Ls)
        U = U1 - U_m * math.exp((y-L/4)/Ls)
    elif y >= L/4 and y < L/2:
        rho = rho2 + rho_m * math.exp((-y+L/4)/Ls)
        U = U2 + U_m * math.exp((-y+L/4)/Ls)
    elif y >= L/2 and y < 3*L/4:
        rho = rho2 + rho_m * math.exp(-(3*L/4-y)/Ls)
        U = U2 + U_m * math.exp(-(3*L/4-y)/Ls)
    else:
        rho = rho1 - rho_m * math.exp(-(y-3*L/4)/Ls)
        U = U1 - U_m * math.exp(-(y-3*L/4)/Ls)
    #
    T = p_inf/(Rgas*rho)
    vy_perturb = 1.0 * math.sin(4*x*math.pi)
    return FlowState(p=p_inf, velx=U, vely=vy_perturb, T=T)
#
blk0 = FluidBlock(grid=grd0, initialState=initial_flow,
                  bcs={'iminus':ExchangeBC(),'iplus':ExchangeBC(),
                       'jminus':ExchangeBC(),'jplus':ExchangeBC()})
#
ref_time = L/abs(U1-U2)
config.max_time = 2.5 * ref_time
config.max_step = 1000000
add_dt_plot(0.0, config.max_time/25)
