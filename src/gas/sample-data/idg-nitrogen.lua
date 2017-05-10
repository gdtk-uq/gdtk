model = "IdealDissociatingGas"

IdealDissociatingGas = {
   molecule = "N2",
   atom = "N",
   W = 28.0, -- g/mole
   T_d = 113200.0, -- K
   rho_d = 130.0, -- g/cm^3
   -- Rate constants follow.
   C1 = 8.5e25, n1 = -2.5,
   C2 = 2.3e29, n2 = -3.5
}
