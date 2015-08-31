config = {
  errTol = 1.000000e-03,
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='rkf'},
}

reaction = {}
reaction[1] = {
  equation = "H2 + I2 <=> 2 HI",
  type = "elementary",
  frc = {model='Arrhenius', A=1.930000000000e+08, n=0.000000, C=2.062000000000e+04 },
  brc = {model='Arrhenius', A=5.470705805110e-07, n=0.000000, C=0.000000000000e+00 },
  ec = {},
  reacIdx = { 0, 1,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 2,},
  prodCoeffs = { 2,},
  anonymousCollider = false,
}
