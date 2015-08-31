config = {
  errTol = 1.000000e-03,
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='rkf'},
}

reaction = {}
reaction[1] = {
  equation = "N2 + N2 <=> N + N + N2",
  type = "elementary",
  frc = {model='Arrhenius', A=7.000000000000e+15, n=-1.600000, C=1.132000000000e+05 },
  brc = {model='Arrhenius', A=1.090000000000e+04, n=-0.500000, C=0.000000000000e+00 },
  ec = {},
  reacIdx = { 0,},
  reacCoeffs = { 2,},
  prodIdx = { 0, 1,},
  prodCoeffs = { 1, 2,},
  anonymousCollider = false,
}
reaction[2] = {
  equation = "N2 + N <=> N + N + N",
  type = "elementary",
  frc = {model='Arrhenius', A=3.000000000000e+16, n=-1.600000, C=1.132000000000e+05 },
  brc = {model='Arrhenius', A=2.320000000000e+09, n=-1.500000, C=0.000000000000e+00 },
  ec = {},
  reacIdx = { 0, 1,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1,},
  prodCoeffs = { 3,},
  anonymousCollider = false,
}
