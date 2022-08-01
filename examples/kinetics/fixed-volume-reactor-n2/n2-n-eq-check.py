# n2-n-eq-check.py
# Usage: python3 n2-n-eq-check.py
from gdtk.gas import GasModel, GasState
gmodel = GasModel('cea-n2-gas-model.lua')
gs = GasState(gmodel)
gs.p = 1.455e+05
gs.T = 6177.424
gs.update_thermo_from_pT()
print("eq gas state=", gs)
print("ceaSavedData=", gs.ceaSavedData)
