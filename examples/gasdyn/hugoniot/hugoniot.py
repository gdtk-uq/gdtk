# hugoniot.py
# Generate the curve in pressure-volume space for a shock jump.
#
# PJ, 2025-02-06,07
#
from gdtk.gas import GasModel, GasState, GasFlow
import math
import numpy as np
import sys

def generate_hugoniot_data(gas_model_file):
    gmodel = GasModel(gas_model_file)
    state1 = GasState(gmodel)
    state1.p = 10.0e3 # Pa
    state1.T = 300.0 # K
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    print("state1: %s" % state1)
    #
    state2 = GasState(gmodel)
    flow = GasFlow(gmodel)
    print("normal shock, given shock speed")
    vss = np.linspace(550.0, 5500.0, 30)
    p_hat = [1.0,]
    v_hat = [1.0,]
    for vs in vss:
        print("vs=%g" % vs)
        v2, vg = flow.normal_shock(state1, vs, state2)
        print("v2=%g vg=%g" % (v2, vg))
        print("state2: %s" % state2)
        p_hat.append(state2.p/state1.p)
        v_hat.append(state1.rho/state2.rho)
    return np.array(v_hat), np.array(p_hat), vss

v_hat_co2, p_hat_co2, vss = generate_hugoniot_data('co2-5sp-eq.lua')
# 'cea-co2-gas-model.lua'
# 'co2-5sp-eq.lua'
v_hat_air, p_hat_air, vss = generate_hugoniot_data('air-5sp-eq.lua')
print("CO2 v_hat=", v_hat_co2)
print("Air v_hat=", v_hat_air)
print("vs=", 0.5*(vss[1:]+vss[:-1]))
print("CO2 delta_v_hat=", (v_hat_co2[1:]-v_hat_co2[:-1]))
print("Air delta_v_hat=", (v_hat_air[1:]-v_hat_air[:-1]))

import matplotlib.pyplot as plt
fig, ax0 = plt.subplots(1,1)
ax0.plot(v_hat_co2, p_hat_co2, label='CO2')
ax0.plot(v_hat_air, p_hat_air, label='air')
ax0.set_xlabel('v2/v1'); ax0.set_ylabel('p2/p1')
# ax0.set_xlim(0,1.0); ax0.set_ylim(0,120)
ax0.legend(loc='upper right')
ax0.set_title('Rankine-Hugoniot shock jump')
if len(sys.argv) == 2 and sys.argv[1] == "--save-figure":
    fig.savefig("hugoniot-curves.png", transparent=True, dpi=300)
else:
    plt.show()
