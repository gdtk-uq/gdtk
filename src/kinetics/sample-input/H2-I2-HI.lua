species = {'H2', 'I2', 'HI'}

H2 = {}
H2.M = 0.00201588
H2.ceaViscosity = {
   nsegments = 1,
   segment0 = {
      A = 0.74553182,
      C = -3257.934,
      B = 43.555109,
      T_upper = 1000,
      T_lower = 200,
      D = 0.13556243,
    }
}
H2.ceaThermCond = {
   nsegments = 1,
   segment0 = {
      A = 1.0059461,
      C = -29792.018,
      B = 279.51262,
      T_upper = 1000,
      T_lower = 200,
      D = 1.1996252,
    },
}
H2.ceaThermoCoeffs = {
   nsegments = 1,
   segment0 = {
    T_upper = 1000,
    T_lower = 200,
    coeffs = {
      40783.2321,
      -800.918604,
      8.21470201,
      -0.01269714457,
      1.753605076e-05,
      -1.20286027e-08,
      3.36809349e-12,
      2682.484665,
      -30.43788844,
    },
  },
}

I2 = {}
I2.M =  0.25380894
I2.ceaViscosity = {
   nsegments = 1,
   segment0 = {
      A = 0.54929498,
      B = -0.36186119e03,
      C = 0.33655931e05,
      D = 0.26154108e01,
      T_upper = 1000,
      T_lower = 200,
    }
}
I2.ceaThermCond = {
   nsegments = 1,
   segment0 = {
      A = 0.29817264,
      B = -0.62470054e03,
      C = 0.33655931e05,
      D = 0.30234067e01,
      T_upper = 1000,
      T_lower = 200,
    },
}
I2.ceaThermoCoeffs = {
   nsegments = 1,
   segment0 = {
    T_upper = 1000,
    T_lower = 200,
    coeffs = {
      -5087.96877,
      -12.4958521,
      4.50421909,
      0.0001370962533,
      -1.390523014e-07,
      1.174813853e-10,
      -2.337541043e-14,
      6213.46981,
      5.58383694,
    },
  },
}

HI = {}
HI.M =  0.12791241
HI.ceaViscosity = {
   nsegments = 1,
   segment0 = {
      A = 0.53718504,
      B = -0.22504609e03,
      C = 0.12416876e05,
      D = 0.27888146e01,
      T_upper = 1000,
      T_lower = 200,
    }
}
HI.ceaThermCond = {
   nsegments = 1,
   segment0 = {
      A = 0.83653272, 
      B = -0.10434645e03,
      C = 0.90075923e04,
      D = -0.38982280,
      T_upper = 1000,
      T_lower = 200,
    },
}
HI.ceaThermoCoeffs = {
   nsegments = 1,
   segment0 = {
    T_upper = 1000,
    T_lower = 200,
    coeffs = {
       18728.8173,
      -343.178884,
      5.95671243,
      -0.0085434396,
      1.454780274e-05,
      -1.049104164e-08,
      2.839734003e-12,
      3682.95072,
      -8.14975609,
    },
  },
}
