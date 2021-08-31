model = "IdealGas"

IdealGas = {
  speciesName = 'air',
  mMass = 0.02896,
  gamma = 1.4,
  entropyRefValues = { 
     s1 = 0.0,
     T1 = 298.15,
     p1 = 101.325e3
  },
  viscosity = {
     model = 'Sutherland',
     mu_ref = 1.716e-5, 
     T_ref = 273.0,
     S = 111.0, 
  },
  thermCondModel = {
     model = 'Sutherland',
     T_ref = 273.0, 
     k_ref = 0.0241, 
     S = 194.0
  }
}
