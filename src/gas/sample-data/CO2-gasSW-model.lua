model = "CO2GasSW"

CO2GasSW= {
  speciesName = 'CO2',
  mMass = 0.04401121121333065,
  gamma = 1.4,
  entropyRefValues = { 
     s1 = 0.0,
     T1 = 298.15,
     p1 = 101.325e3
  },
  sutherlandVisc = {
     mu_ref = 1.716e-5, 
     T_ref = 273.0,
     S = 111.0, 
  },
  sutherlandThermCond = {
     T_ref = 273.0, 
     k_ref = 0.0241, 
     S = 194.0
  },
  LUTfilenames = {
    p_rhoe_file = 'LUT/P_rhoe_Tree.dat',
    a_rhoe_file = 'LUT/a_rhoe_Tree.dat',
    T_rhoe_file = 'LUT/T_rhoe_Tree.dat',
  },
}
