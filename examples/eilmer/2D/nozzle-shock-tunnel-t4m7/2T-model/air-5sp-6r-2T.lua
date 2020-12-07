species = {[0]='N2', [1]='O2', [2]='N', [3]='O', [4]='NO', }
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
  frc = {model='Park', A=3.610000000000e+12, n=-1.000000, C=5.940000000000e+04, s=0.500000 },
  brc = {model='Arrhenius', A=3.010000000000e+03, n=-0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 1,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 3,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [0]=2.000000e+00,
    [1]=9.000000e+00,
    [2]=1.000000e+00,
    [3]=2.500000e+01,
    [4]=1.000000e+00,
  },
}

reaction[2] = {
  equation = "N2 + M <=> N + N + M",
  type = "anonymous_collider",
  frc = {model='Park', A=1.920000000000e+11, n=-0.500000, C=1.131000000000e+05, s=0.500000 },
  brc = {model='Arrhenius', A=1.090000000000e+04, n=-0.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 2,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [0]=2.500000e+00,
    [1]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
  },
}

reaction[3] = {
  equation = "N2 + N <=> N + N + N",
  type = "elementary",
  frc = {model='Park', A=4.150000000000e+16, n=-1.500000, C=1.131000000000e+05, s=0.500000 },
  brc = {model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2,},
  prodCoeffs = { 3.000000e+00,},
}

reaction[4] = {
  equation = "NO + M <=> N + O + M",
  type = "anonymous_collider",
  frc = {model='Park', A=3.970000000000e+14, n=-1.500000, C=7.560000000000e+04, s=0.500000 },
  brc = {model='Arrhenius', A=1.010000000000e+08, n=-1.500000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 4,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 2, 3,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.000000e+00,
    [2]=2.000000e+01,
    [3]=2.000000e+01,
    [4]=2.000000e+01,
  },
}

reaction[5] = {
  equation = "NO + O <=> O2 + N",
  type = "elementary",
  frc = {model='Arrhenius', A=3.180000000000e+03, n=1.000000, C=1.970000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=9.630000000000e+05, n=0.500000, C=3.600000000000e+03, rctIndex=-1},
  ec = {},
  reacIdx = { 3, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1, 2,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[6] = {
  equation = "N2 + O <=> NO + N",
  type = "elementary",
  frc = {model='Arrhenius', A=6.750000000000e+07, n=0.000000, C=3.750000000000e+04, rctIndex=-1},
  brc = {model='Arrhenius', A=1.500000000000e+07, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  ec = {},
  reacIdx = { 0, 3,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

