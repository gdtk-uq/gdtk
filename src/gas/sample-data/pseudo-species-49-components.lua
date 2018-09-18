-- Author: Rowan J. Gollan
-- Date: 2018-08-06
--
-- Demo file for defining pseudo-species components in a 
-- pseudo-species gas model.

model = "PseudoSpeciesGas"
number_pseudo_species = 49
pseudo_species = {}
number_parents_species = 2
parents = {"N","N2",}
        
pseudo_species[0] = {
   name = "N_e2p3Minus4S",
   type = "single_state",
   M = 0.014,
   DOF_base_mode = 3,
   energy = 0.0,
   parentIdx = 0
}
    
pseudo_species[1] = {
   name = "N2_eX1SIGGPlus_v0",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -9.754153395531866,
   parentIdx = 1
}
    
pseudo_species[2] = {
   name = "N2_eX1SIGGPlus_v1",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -9.465432276631866,
   parentIdx = 1
}
    
pseudo_species[3] = {
   name = "N2_eX1SIGGPlus_v2",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -9.180261285931866,
   parentIdx = 1
}
    
pseudo_species[4] = {
   name = "N2_eX1SIGGPlus_v3",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -8.898640423431866,
   parentIdx = 1
}
    
pseudo_species[5] = {
   name = "N2_eX1SIGGPlus_v4",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -8.620569689031866,
   parentIdx = 1
}
    
pseudo_species[6] = {
   name = "N2_eX1SIGGPlus_v5",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -8.346049082931867,
   parentIdx = 1
}
    
pseudo_species[7] = {
   name = "N2_eX1SIGGPlus_v6",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -8.075078604931866,
   parentIdx = 1
}
    
pseudo_species[8] = {
   name = "N2_eX1SIGGPlus_v7",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -7.807658255131866,
   parentIdx = 1
}
    
pseudo_species[9] = {
   name = "N2_eX1SIGGPlus_v8",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -7.543788033531866,
   parentIdx = 1
}
    
pseudo_species[10] = {
   name = "N2_eX1SIGGPlus_v9",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -7.283467940031866,
   parentIdx = 1
}
    
pseudo_species[11] = {
   name = "N2_eX1SIGGPlus_v10",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -7.026697974831865,
   parentIdx = 1
}
    
pseudo_species[12] = {
   name = "N2_eX1SIGGPlus_v11",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -6.773478137731866,
   parentIdx = 1
}
    
pseudo_species[13] = {
   name = "N2_eX1SIGGPlus_v12",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -6.523808428831866,
   parentIdx = 1
}
    
pseudo_species[14] = {
   name = "N2_eX1SIGGPlus_v13",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -6.277688848131866,
   parentIdx = 1
}
    
pseudo_species[15] = {
   name = "N2_eX1SIGGPlus_v14",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -6.035119395631865,
   parentIdx = 1
}
    
pseudo_species[16] = {
   name = "N2_eX1SIGGPlus_v15",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -5.796100071331866,
   parentIdx = 1
}
    
pseudo_species[17] = {
   name = "N2_eX1SIGGPlus_v16",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -5.560630875131866,
   parentIdx = 1
}
    
pseudo_species[18] = {
   name = "N2_eX1SIGGPlus_v17",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -5.3287118071318655,
   parentIdx = 1
}
    
pseudo_species[19] = {
   name = "N2_eX1SIGGPlus_v18",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -5.100342867331865,
   parentIdx = 1
}
    
pseudo_species[20] = {
   name = "N2_eX1SIGGPlus_v19",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -4.875524055731866,
   parentIdx = 1
}
    
pseudo_species[21] = {
   name = "N2_eX1SIGGPlus_v20",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -4.654255372331866,
   parentIdx = 1
}
    
pseudo_species[22] = {
   name = "N2_eX1SIGGPlus_v21",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -4.4365368171318655,
   parentIdx = 1
}
    
pseudo_species[23] = {
   name = "N2_eX1SIGGPlus_v22",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -4.222368390031866,
   parentIdx = 1
}
    
pseudo_species[24] = {
   name = "N2_eX1SIGGPlus_v23",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -4.011750091231866,
   parentIdx = 1
}
    
pseudo_species[25] = {
   name = "N2_eX1SIGGPlus_v24",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -3.8046819205318654,
   parentIdx = 1
}
    
pseudo_species[26] = {
   name = "N2_eX1SIGGPlus_v25",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -3.6011638780318656,
   parentIdx = 1
}
    
pseudo_species[27] = {
   name = "N2_eX1SIGGPlus_v26",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -3.4011959636318663,
   parentIdx = 1
}
    
pseudo_species[28] = {
   name = "N2_eX1SIGGPlus_v27",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -3.204778177531866,
   parentIdx = 1
}
    
pseudo_species[29] = {
   name = "N2_eX1SIGGPlus_v28",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -3.011910519531865,
   parentIdx = 1
}
    
pseudo_species[30] = {
   name = "N2_eX1SIGGPlus_v29",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -2.8225929898318665,
   parentIdx = 1
}
    
pseudo_species[31] = {
   name = "N2_eX1SIGGPlus_v30",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -2.636825588231865,
   parentIdx = 1
}
    
pseudo_species[32] = {
   name = "N2_eX1SIGGPlus_v31",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -2.454608314831866,
   parentIdx = 1
}
    
pseudo_species[33] = {
   name = "N2_eX1SIGGPlus_v32",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -2.2759411695318654,
   parentIdx = 1
}
    
pseudo_species[34] = {
   name = "N2_eX1SIGGPlus_v33",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -2.1008241525318656,
   parentIdx = 1
}
    
pseudo_species[35] = {
   name = "N2_eX1SIGGPlus_v34",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -1.9292572636318663,
   parentIdx = 1
}
    
pseudo_species[36] = {
   name = "N2_eX1SIGGPlus_v35",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -1.7612405029318658,
   parentIdx = 1
}
    
pseudo_species[37] = {
   name = "N2_eX1SIGGPlus_v36",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -1.596773870531866,
   parentIdx = 1
}
    
pseudo_species[38] = {
   name = "N2_eX1SIGGPlus_v37",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -1.4358573661318665,
   parentIdx = 1
}
    
pseudo_species[39] = {
   name = "N2_eX1SIGGPlus_v38",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -1.278490990031866,
   parentIdx = 1
}
    
pseudo_species[40] = {
   name = "N2_eX1SIGGPlus_v39",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -1.124674742131866,
   parentIdx = 1
}
    
pseudo_species[41] = {
   name = "N2_eX1SIGGPlus_v40",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.9744086223318664,
   parentIdx = 1
}
    
pseudo_species[42] = {
   name = "N2_eX1SIGGPlus_v41",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.8276926307318657,
   parentIdx = 1
}
    
pseudo_species[43] = {
   name = "N2_eX1SIGGPlus_v42",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.6845267673318656,
   parentIdx = 1
}
    
pseudo_species[44] = {
   name = "N2_eX1SIGGPlus_v43",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.544911032131866,
   parentIdx = 1
}
    
pseudo_species[45] = {
   name = "N2_eX1SIGGPlus_v44",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.40884542513186517,
   parentIdx = 1
}
    
pseudo_species[46] = {
   name = "N2_eX1SIGGPlus_v45",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.2763299462318667,
   parentIdx = 1
}
    
pseudo_species[47] = {
   name = "N2_eX1SIGGPlus_v46",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.14736459563186521,
   parentIdx = 1
}
    
pseudo_species[48] = {
   name = "N2_eX1SIGGPlus_v47",
   type = "single_state",
   M = 0.028014,
   DOF_base_mode = 5,
   energy = -0.021949373131866068,
   parentIdx = 1
}

---------- Transport
db = {}
db['N2'] = {}
db['N2'].atomicConstituents = { N=2, }
db['N2'].charge = 0
db['N2'].M = 0.02801340
db['N2'].sigma = 3.62100000
db['N2'].epsilon = 97.53000000


db['N2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.2526577e-01,
      B = -3.1779652e+01,
      C = -1.6407983e+03,
      D =  1.7454992e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.7395209e-01,
      B =  5.6152222e+02,
      C = -1.7394809e+05,
      D = -3.9335958e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.8503551e-01,
      B =  9.0902171e+02,
      C = -7.3129061e+05,
      D = -5.3503838e-01,
   },
}
db['N2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  8.5439436e-01,
      B =  1.0573224e+02,
      C = -1.2347848e+04,
      D =  4.7793128e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.8407146e-01,
      B =  1.3357293e+02,
      C = -1.1429640e+04,
      D =  2.4417019e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  2.4176185e+00,
      B =  8.0477749e+03,
      C =  3.1055802e+06,
      D = -1.4517761e+01,
   },
}


db['N'] = {}
db['N'].atomicConstituents = { N=1, }
db['N'].charge = 0
db['N'].M = 0.01400670
db['N'].sigma = 3.29800000
db['N'].epsilon = 71.40000000
db['N'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 200.0,
    T_upper = 1000.0,
    coeffs = {
       0.000000000e+00,
       0.000000000e+00,
       2.500000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       5.610463780e+04,
       4.193905036e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       8.876501380e+04,
      -1.071231500e+02,
       2.362188287e+00,
       2.916720081e-04,
      -1.729515100e-07,
       4.012657880e-11,
      -2.677227571e-15,
       5.697351330e+04,
       4.865231506e+00,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       5.475181050e+08,
      -3.107574980e+05,
       6.916782740e+01,
      -6.847988130e-03,
       3.827572400e-07,
      -1.098367709e-11,
       1.277986024e-16,
       2.550585618e+06,
      -5.848769753e+02,
    }
  },
}
db['N'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3724737e-01,
      B =  4.3997150e+02,
      C = -1.7450753e+05,
      D =  1.0365689e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.9986588e-01,
      B =  1.4112801e+03,
      C = -1.8200478e+06,
      D = -5.5811716e-01,
   },
}
db['N'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3771661e-01,
      B =  4.4243270e+02,
      C = -1.7578446e+05,
      D =  8.9942915e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  9.0001710e-01,
      B =  1.4141175e+03,
      C = -1.8262403e+06,
      D =  2.4048513e-01,
   },
}

