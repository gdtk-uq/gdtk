# n2-n-eq-check.py
# Usage: python3 n2-n-eq-check.py
from gdtk.gas import GasModel, GasState
gmodel = GasModel('cea-n2-gas-model.lua')
gs = GasState(gmodel)
gs.p = 1.0e+05
gs.T = 6001.864
gs.update_thermo_from_pT()
print("eq gas state=", gs)
print("ceaSavedData=", gs.ceaSavedData)
