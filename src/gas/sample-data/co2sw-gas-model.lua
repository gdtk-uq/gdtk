model = 'CO2GasSW'
CO2GasSW = {
   speciesName = 'CO2',
   mMass = 0.04401000,
   gamma = 1.30000000,
   entropyRefValues = {
      s1 = 0.00000000e+00,
      T1 = 298.15000000,
      p1 = 1.01325000e+05,
   },
   sutherlandVisc = {
      mu_ref = 1.48000000e-05,
      T_ref = 293.15000000,
      S = 240.00000000,
   },
   sutherlandThermCond = {
      k_ref = 2.41000000e-02,
      T_ref = 273.00000000,
      S = 194.00000000,
   },
   LUTfilenames = {
      p_rhoe_file = './LUT/P_rhoe_Tree.dat',
      a_rhoe_file = './LUT/a_rhoe_Tree.dat',
      T_rhoe_file = './LUT/T_rhoe_Tree.dat',
      e_rho_sat_file = './LUT/e_rho_sat_table.dat',
      rho_sh_file = './LUT/rho_sh_Tree.dat',
      T_sh_file = './LUT/T_sh_Tree.dat',
      lookup_hsFlag = 1,
      lookup_rhoeFlag = 1,
   }
}
