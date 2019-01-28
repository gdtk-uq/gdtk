config = {
  tempLimits = {lower=300.000000, upper=30000.000000},
  odeStep = {method='rkf', errTol=1.000000e-03},
  tightTempCoupling = false,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
reaction[1] = {
  equation = "O2 + M <=> O + O + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=3.610000000000e+12, n=-1.000000, C=5.940000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=3.010000000000e+03, n=-0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 3,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 1,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.500000e+01,
    [2]=2.000000e+00,
    [3]=9.000000e+00,
    [4]=1.000000e+00,
  },
}

reaction[2] = {
  equation = "N2 + M <=> N + N + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=1.920000000000e+11, n=-0.500000, C=1.131000000000e+05, rctIndex=-1},
  brc = {model='Arrhenius', A=1.090000000000e+04, n=-0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 2,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 0,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [1]=1.000000e+00,
    [2]=2.500000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
  },
}

reaction[3] = {
  equation = "N2 + N <=> N + N + N",
  type = "elementary",
  frc = {model='Arrhenius', A=4.150000000000e+16, n=-1.500000, C=1.131000000000e+05, rctIndex=-1},
  brc = {model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0,},
  prodCoeffs = { 3.000000e+00,},
}

reaction[4] = {
  equation = "NO + M <=> N + O + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=3.970000000000e+14, n=-1.500000, C=7.560000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.010000000000e+08, n=-1.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 4,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 0, 1,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  efficiencies = {
    [0]=2.000000e+01,
    [1]=2.000000e+01,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=2.000000e+01,
  },
}

reaction[5] = {
  equation = "NO + O <=> O2 + N",
  type = "elementary",
  frc = {model='Arrhenius', A=3.180000000000e+03, n=1.000000, C=1.970000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=9.630000000000e+05, n=0.500000, C=3.600000000000e+03, rctIndex=-1},
  ec = {},
  reacIdx = { 1, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 3,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[6] = {
  equation = "N2 + O <=> NO + N",
  type = "elementary",
  frc = {model='Arrhenius', A=6.750000000000e+07, n=0.000000, C=3.750000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.500000000000e+07, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[7] = {
  equation = "N + O <=> NO+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=9.030000000000e+03, n=0.500000, C=3.240000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.800000000000e+13, n=-1.000000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0, 1,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 9, 10,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[8] = {
  equation = "O + e- <=> O+ + e- + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=3.600000000000e+25, n=-2.910000, C=1.580000000000e+05, rctIndex=-1},
  brc = {model='Arrhenius', A=2.200000000000e+28, n=-4.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1, 10,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 6, 10,},
  prodCoeffs = { 1.000000e+00, 2.000000e+00,},
}

reaction[9] = {
  equation = "N + e- <=> N+ + e- + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=1.100000000000e+26, n=-3.140000, C=1.690000000000e+05, rctIndex=-1},
  brc = {model='Arrhenius', A=2.200000000000e+28, n=-4.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0, 10,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5, 10,},
  prodCoeffs = { 1.000000e+00, 2.000000e+00,},
}

reaction[10] = {
  equation = "O + O <=> O2+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=1.600000000000e+11, n=-0.980000, C=8.080000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=8.020000000000e+15, n=-1.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 8, 10,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[11] = {
  equation = "O + O2+ <=> O2 + O+",
  type = "elementary",
  frc = {model='Arrhenius', A=2.920000000000e+12, n=-1.110000, C=2.800000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=7.800000000000e+05, n=0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1, 8,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 6,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[12] = {
  equation = "N2 + N+ <=> N + N2+",
  type = "elementary",
  frc = {model='Arrhenius', A=2.020000000000e+05, n=0.810000, C=1.300000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=7.800000000000e+05, n=0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 2, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[13] = {
  equation = "N + N <=> N2+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=1.400000000000e+07, n=0.000000, C=6.780000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.500000000000e+16, n=-1.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 7, 10,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[14] = {
  equation = "O2 + N2 <=> NO + NO+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=1.380000000000e+14, n=-1.840000, C=1.410000000000e+05, rctIndex=-1},
  brc = {model='Arrhenius', A=1.000000000000e+12, n=-2.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 2, 3,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 4, 9, 10,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00, 1.000000e+00,},
}

reaction[15] = {
  equation = "NO + M <=> NO+ + e- + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=2.200000000000e+09, n=-0.350000, C=1.008000000000e+05, rctIndex=-1},
  brc = {model='Arrhenius', A=2.200000000000e+14, n=-2.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 4,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 9, 10,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  efficiencies = {
    [2]=1.000000e+00,
    [3]=4.000000e+00,
  },
}

reaction[16] = {
  equation = "O + NO+ <=> NO + O+",
  type = "elementary",
  frc = {model='Arrhenius', A=3.630000000000e+09, n=-0.600000, C=5.080000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.500000000000e+07, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1, 9,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 4, 6,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[17] = {
  equation = "N2 + O+ <=> O + N2+",
  type = "elementary",
  frc = {model='Arrhenius', A=3.400000000000e+13, n=-2.000000, C=2.300000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=2.480000000000e+13, n=-2.200000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 2, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[18] = {
  equation = "N + NO+ <=> NO + N+",
  type = "elementary",
  frc = {model='Arrhenius', A=1.000000000000e+13, n=-0.930000, C=6.100000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=4.800000000000e+08, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0, 9,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 4, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[19] = {
  equation = "O2 + NO+ <=> NO + O2+",
  type = "elementary",
  frc = {model='Arrhenius', A=1.800000000000e+09, n=0.170000, C=3.300000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.800000000000e+07, n=0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 3, 9,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 4, 8,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[20] = {
  equation = "O + NO+ <=> O2 + N+",
  type = "elementary",
  frc = {model='Arrhenius', A=1.340000000000e+07, n=0.310000, C=7.727000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.000000000000e+08, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1, 9,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

