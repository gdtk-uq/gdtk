# sample input for Puffin flow solver.
# This approximates the cone20 test case from Eilmer.
# PJ 2022-01-21

init_gas_model("ideal-air-gas-model.lua")
gas = GasState(config.gmodel)
gas.p = 100.0e3
gas.T = 300.0

config.axisymmetric = True
config.max_step_relax = 20

def lower_y(x):
    return math.tan(math.radians(20.0))*(x-0.2) if x > 0.2 else 0.0
def upper_y(x): return 1.0
def lower_bc(x): return 0
def upper_bc(x): return 0

st1 = StreamTube(gas=gas, velx=1000.0, vely=0.0,
                 y0=lower_y, y1=upper_y, bc0=lower_bc, bc1=upper_bc,
                 ncells=100)
