-- 11 Species Air Reaction File, from:
-- "Modification of chemical-kinetic parameters for 11-air species in re-entry flows"
--  Jae Gang Kim and Sung Min Jo, International Journal of Heat and Mass Transfer, Volume 169, 2021
--
-- Notes:
--    - Note that in the input file (here), the reactions are specified in CGS units,
--      primarily cm. prep-chem converts these into SI units automatically for the 
--      machine readable reaction file that gets fed into the lua interpreter at runtime.
--    - Table 3 does not appear to list the units for the pre-exponential factor,
--      but cross referencing with reference [42] reveals they are in cm^3/mol/s
-- @author: Nick Gibbons (21/04/30)

Config{
   tempLimits = {lower=300.0, upper=30000.0},
   odeStep = {method='rkf', errTol=1.000000e-09},
   tightTempCoupling = true,
}

energyModes = {["vibration"]=0, ["electron-electronic"]=1}

-- Park-type Reactions are evaluated as:
-- T = pow(Q.T, s)*pow(Q.T_modes[0], 1.0 - s);
-- return A*pow(T, n)*exp(-C/T);

-- Dissociation Reactions (other than Electrons) from Table 3
Reaction{'N2 + N2 <=> N + N + N2',
   fr={'Park', A=1.216e+20, n=-1.214, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + O2 <=> N + N + O2',
   fr={'Park', A=7.000e+21, n=-1.600, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + NO <=> N + N + NO',
   fr={'Park', A=7.000e+21, n=-1.600, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + N <=> N + N + N',
   fr={'Park', A=3.591e+20, n=-1.226, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + O <=> N + N + O',
   fr={'Park', A=3.000e+22, n=-1.600, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + N2+ <=> N + N + N2+',
   fr={'Park', A=1.216e+20, n=-1.214, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + O2+ <=> N + N + O2+',
   fr={'Park', A=7.000e+21, n=-1.600, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + NO+ <=> N + N + NO+',
   fr={'Park', A=7.000e+21, n=-1.600, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + N+ <=> N + N + N+',
   fr={'Park', A=3.591e+20, n=-1.226, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'N2 + O+ <=> N + N + O+',
   fr={'Park', A=3.000e+22, n=-1.600, C=113200, s=0.5}, 
   label="N2diss"
}

Reaction{'O2 + N2 <=> O + O + N2',
   fr={'Park', A=3.354e+15, n=-0.2726, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + O2 <=> O + O + O2',
   fr={'Park', A=1.117e+25, n=-2.585, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + NO <=> O + O + NO',
   fr={'Park', A=3.354e+15, n=-0.2726, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + N <=> O + O + N',
   fr={'Park', A=1.000e+22, n=-1.500, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + O <=> O + O + O',
   fr={'Park', A=3.000e+21, n=-1.500, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + N2+ <=> O + O + N2+',
   fr={'Park', A=3.354e+15, n=-0.2726, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + O2+ <=> O + O + O2+',
   fr={'Park', A=1.117e+25, n=-2.585, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + NO+ <=> O + O + NO+',
   fr={'Park', A=3.354e+15, n=-0.2726, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + N+ <=> O + O + N+',
   fr={'Park', A=1.000e+22, n=-1.500, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'O2 + O+ <=> O + O + O+',
   fr={'Park', A=3.000e+21, n=-1.500, C=59500, s=0.5}, 
   label="O2diss"
}

Reaction{'NO + M <=> N + O + M',
   fr={'Park', A=1.450e+15, n=0.0, C=75200.0, s=0.5},
   efficiencies={['N2']=1.0, ['N2+']=1.0,
                 ['O2']=1.0, ['O2+']=1.0,
                 ['NO']=0.664827586,['NO+']=0.664827586,
                 ['N'] =0.664827586, ['N+']=0.664827586,
                 ['O'] =0.664827586, ['O+']=0.664827586,
                 ['e-']=0.0},
   label="NOdiss"
}

-- The rest of the mechanism is the same as Park 1993
-- Dissociation Reactions (electron impact)
Reaction{ 'N2 + e- <=> N + N + e-',
   fr={'Park', A=1.2e25, n=-1.60, C=113200.0, s=0.5, rateControllingTemperature="electron-electronic"},
   br = {"fromEqConst", rateControllingTemperature="electron-electronic"},
   label="N2EIdiss"
}
-- For some reason the efficients of O2 + e- and NO + e- are zero ?

-- Arrhenius reaction rates are evaluated using the translation temperature only:
-- return A*pow(T, n)*exp(-C/T);

-- NO Exchange Reactions
Reaction{ 'NO + O <=> O2 + N',
   fr={'Arrhenius', A=8.4e12, n=0.0, C=19450.0},
   label="exchange"
}

Reaction{'N2 + O <=> NO + N',
   fr={'Arrhenius', A=6.4e17, n=-1.0, C=38400.0},
   label="exchange"
}

-- Associative ionization reactions
Reaction{'N + O <=> NO+ + e-',
   fr={'Arrhenius', A=8.8e8, n=1.0, C=31900.0},
   br={'fromEqConst', rateControllingTemperature="electron-electronic"},
   label="ass_ion"
}

Reaction{'O + O <=> O2+ + e-',
   fr={'Arrhenius', A=7.1e2, n=2.7, C=80600.0},
   br={'fromEqConst', rateControllingTemperature='electron-electronic'},
   label="ass_ion"
}

Reaction{'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.4e7, n=1.5, C=67500.0},
   br={'fromEqConst', rateControllingTemperature='electron-electronic'},
   label="ass_ion"
}

-- Charge exchange reactions
Reaction{'NO+ + O <=> N+ + O2',
   fr={'Arrhenius', A=1.0e12, n=0.5, C=77200.0},
   label="charge_exchange"
}

Reaction{'N+ + N2 <=> N2+ + N',
   fr={'Arrhenius', A=1.0e12, n=0.5, C=12200.0},
   label="charge_exchange"
}

Reaction{'O2+ + N <=> N+ + O2',
   fr={'Arrhenius', A=8.7e13, n=0.14, C=28600.0},
   label="charge_exchange"
}

Reaction{'O+ + NO <=> N+ + O2',
   fr={'Arrhenius', A=1.4e5, n=1.9, C=26600.0},
   label="charge_exchange"
}

Reaction{'O2+ + N2 <=> N2+ + O2',
   fr={'Arrhenius', A=9.9e12, n=0.0, C=40700.0},
   label="charge_exchange"
}

Reaction{'O2+ + O <=> O+ + O2',
   fr={'Arrhenius', A=4.0e12, n=-0.09, C=18000.0},
   label="charge_exchange"
}

Reaction{'NO+ + N <=> O+ + N2',
   fr={'Arrhenius', A=3.4e13, n=-1.08, C=12800.0},
   label="charge_exchange"
}

Reaction{'NO+ + O2 <=> O2+ + NO',
   fr={'Arrhenius', A=2.4e13,  n=0.41, C=32600.0},
   label="charge_exchange"
}

Reaction{'NO+ + O <=> O2+ + N',
   fr={'Arrhenius', A=7.2e12,  n=0.29, C=48600.0},
   label="charge_exchange"
}

Reaction{'O+ + N2 <=> N2+ + O',
   fr={'Arrhenius', A=9.1e11,  n=0.36, C=22800.0},
   label="charge_exchange"
}

Reaction{'NO+ + N <=> N2+ + O',
   fr={'Arrhenius', A=7.2e13,  n=0.00, C=35500.0},
   label="charge_exchange"
}

-- Electron-impact ionization reactions
Reaction{'O + e- <=> O+ + e- + e-',
   fr={'Park', A=3.9e33, n=-3.78, C=158500.0, s=0.0, rateControllingTemperature='electron-electronic'},
   br={'fromEqConst', rateControllingTemperature='electron-electronic'},
   label="EII"
}

Reaction{'N + e- <=> N+ + e- + e-',
   fr={'Park', A=2.5e34, n=-3.82, C=168600.0, s=0.0, rateControllingTemperature='electron-electronic'},
   br={'fromEqConst', rateControllingTemperature='electron-electronic'},
   label="EII"
}

-- Radiative recombination reactions
-- Dan had these commented out for some reason. For the sake of comparison we'll leave them out too.
-- Reaction{'O+ + e- <=> O',
--     fr={'Park', A=1.07e11, n=-0.52, C=0.0, s=0.0},
-- }

-- Reaction{'N+ + e- <=> N',
--     fr={'Park', A=1.52e11, n=-0.48, C=0.0, s=0.0},
-- }
