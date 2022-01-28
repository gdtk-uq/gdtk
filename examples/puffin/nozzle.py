# Simple conical nozzle to try axisymmetric terms.
# PJ 2022-01-29

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
gas1.p = 100.0e3
gas1.T = 300.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 1.01
V1 = M1 * gas1.a
print("V1=", V1)

config.axisymmetric = True
config.max_step_relax = 80
config.max_x = 1.0

def lower_y(x): return 0.0
def upper_y(x): return 0.1 + math.tan(math.radians(5.0))*x
def lower_bc(x): return 0
def upper_bc(x): return 0

y_throat = upper_y(0.0)
y_exit = upper_y(config.max_x)
area_ratio = (y_exit/y_throat)**2
print("y_throat=", y_throat, "y_exit=", y_exit, "area_ratio=", area_ratio)

st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=lower_y, y1=upper_y, bc0=lower_bc, bc1=upper_bc,
                 ncells=40)
