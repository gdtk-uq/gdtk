# sample input for Puffin flow solver.
# If axisymmetric:
#     This approximates a Mach 2.0 flow over a cone,
#     with 10-degree semi-vertex angle.
#     From NACA Report 1135, we get the approximate values.
#     Chart 4 post-shock Mach number 1.64
#     Chart 5 shock angle 31.3 degrees
#     Chart 6 (pc-p1)/q1 = 0.105
#
# PJ 2022-01-21

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
gas1.p = 100.0e3
gas1.T = 300.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 2.0
V1 = M1 * gas1.a
print("V1=", V1)
q1 = 0.5*gas1.rho*V1*V1
print("q1=", q1)
pcone = gas1.p + 0.105*q1
print("Expected pcone=", pcone)

# config.axisymmetric = True
config.max_step_relax = 20

def lower_y(x):
    return math.tan(math.radians(10.0))*(x-0.2) if x > 0.2 else 0.0
def upper_y(x): return 1.0
def lower_bc(x): return 0
def upper_bc(x): return 0

st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=lower_y, y1=upper_y, bc0=lower_bc, bc1=upper_bc,
                 ncells=100)
