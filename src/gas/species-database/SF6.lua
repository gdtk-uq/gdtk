db.SF6 = {}--NONE OF SF6 DATA HAS BEEN UDPATED THIS IS MERELY A PLACED HOLDER
db.SF6.atomicConstituents = {}
db.SF6.charge = 0
db.SF6.M = {
   value = 0.146055,
   units = 'kg/mol',
}
db.SF6.gamma = {
   value = 1.3,
   note = "valid at low temperatures"
}
db.SF6.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.SF6.sutherlandVisc = {
   mu_ref = 14.8e-6, 
   T_ref = 293.15,
   S = 240.0,
   reference = "Crane Company (1988) - Flow of fluids through valves, fittings and pipes"
}
db.SF6.sutherlandThermCond = {
   T_ref = 273.0, --these have not been updated
   k_ref = 0.0241, --these have not been updated
   S = 194.0,--these have not been updated
   reference = "Table 1-3, White (2006)"
}

db.SF6.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      3.309526740e+05,
     -4.737685050e+03,
      2.247738068e+01,
      1.046954309e-02,
     -2.560641961e-05,
      2.153716967e-08,
     -6.516098960e-12,
     -1.255360583e+05,
     -1.091760145e+02
   },
   segment1 = {
     -7.306726500e+05,
     -6.367056550e+02,
      1.947442853e+01,
     -1.894325671e-04,
      4.178722830e-08,
     -4.783744950e-12,
      2.213516129e-16,
     -1.510609837e+05,
     -8.147574587e+01
   }
}
db.SF6.Hf = {
   value = -1219400.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
