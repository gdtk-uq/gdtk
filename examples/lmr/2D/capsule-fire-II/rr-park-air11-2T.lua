-- 11 Species Air Reaction File, from:
-- "Review of Chemical Kinetic Problems of Future Nasa Missions, I: Earth Entries"
-- Chul Park, Journal of Thermophysics and Heat Transfer, Vol. 7, No. 3 (1993)
-- Notes:
--    - Note that in the input file (here), the reactions are specified in CGS units,
--      primarily cm. prep-chem converts these into SI units automatically for the 
--      machine readable reaction file that gets fed into the lua interpreter at runtime.
--    - This file was ported from an eilmer 3 file written by Daniel Potter
-- @author: Nick Gibbons

Config{
   tempLimits = {lower=300.0, upper=30000.0},
   odeStep = {method='rkf', errTol=1.000000e-09},
   tightTempCoupling = true,
}

-- Park Reactions are evaluated as:
-- T = pow(Q.T, s)*pow(Q.T_modes[0], 1.0 - s);
-- return A*pow(T, n)*exp(-C/T);

-- Dissociation Reactions
Reaction{ 'N2 + M <=> N + N + M',
   fr={'Park', A=7.0e21, n=-1.60, C=113200.0, s=0.5},
   efficiencies={['N2']=1.0,   ['N2+']=1.0,
                 ['O2']=1.0,   ['O2+']=1.0,
                 ['NO']=1.0,   ['NO+']=1.0,
                 ['N']=4.2857, ['N+']=4.2857,
                 ['O']=4.2857, ['O+']=4.2857,
                 ['e-']=1714.2857142857142},
}

Reaction{'O2 + M <=> O + O + M',
   fr={'Park', A=2.0e21, n=-1.5, C=59500.0, s=0.5},
   efficiencies={['N2']=1.0, ['N2+']=1.0,
                 ['O2']=1.0, ['O2+']=1.0,
                 ['NO']=1.0, ['NO+']=1.0,
                 ['N']=5.0,  ['N+']=5.0,
                 ['O']=5.0,  ['O+']=5.0,
                 ['e-']=0.0},
}

-- This reaction has a typo in the 1993 paper, I assume it should be NO+M -> N+O+M
-- Also for some reason Dan Potter has slightly different efficiencies for the atomic species.
-- I've left these in for the sake of comparison but we should think about putting them back
Reaction{'NO + M <=> N + O + M',
   fr={'Park', A=5.0e15, n=0.0, C=75500.0, s=0.5},
   efficiencies={['N2']=1.0, ['N2+']=1.0,
                 ['O2']=1.0, ['O2+']=1.0,
                 ['NO']=1.0, ['NO+']=1.0,
                 ['N']=20.0, ['N+']=20.0,
                 ['O']=20.0, ['O+']=20.0,
                 ['e-']=0.0},
}


-- Arrhenius reaction rates are evaluated using the translation temperature only:
-- return A*pow(T, n)*exp(-C/T);

-- NO Exchange Reactions
Reaction{ 'NO + O <=> O2 + N',
   fr={'Arrhenius', A=8.4e12, n=0.0, C=19450.0},
}

Reaction{'N2 + O <=> NO + N',
   fr={'Arrhenius', A=6.4e17, n=-1.0, C=38400.0},
}

-- Associative ionization reactions
Reaction{'N + O <=> NO+ + e-',
   fr={'Arrhenius', A=8.8e8, n=1.0, C=31900.0},
}

Reaction{'O + O <=> O2+ + e-',
   fr={'Arrhenius', A=7.1e2, n=2.7, C=80600.0},
}

Reaction{'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.4e7, n=1.5, C=67500.0},
}

-- Charge exchange reactions
Reaction{'NO+ + O <=> N+ + O2',
   fr={'Arrhenius', A=1.0e12, n=0.5, C=77200.0},
}

Reaction{'N+ + N2 <=> N2+ + N',
   fr={'Arrhenius', A=1.0e12, n=0.5, C=12200.0},
}

Reaction{'O2+ + N <=> N+ + O2',
   fr={'Arrhenius', A=8.7e13, n=0.14, C=28600.0},
}

Reaction{'O+ + NO <=> N+ + O2',
   fr={'Arrhenius', A=1.4e5, n=1.9, C=26600.0},
}

Reaction{'O2+ + N2 <=> N2+ + O2',
   fr={'Arrhenius', A=9.9e12, n=0.0, C=40700.0},
}

Reaction{'O2+ + O <=> O+ + O2',
   fr={'Arrhenius', A=4.0e12, n=-0.09, C=18000.0},
}

Reaction{'NO+ + N <=> O+ + N2',
   fr={'Arrhenius', A=3.4e13, n=-1.08, C=12800.0},
}

Reaction{'NO+ + O2 <=> O2+ + NO',
   fr={'Arrhenius', A=2.4e13,  n=0.41, C=32600.0},
}

Reaction{'NO+ + O <=> O2+ + N',
   fr={'Arrhenius', A=7.2e12,  n=0.29, C=48600.0},
}

Reaction{'O+ + N2 <=> N2+ + O',
   fr={'Arrhenius', A=9.1e11,  n=0.36, C=22800.0},
}

Reaction{'NO+ + N <=> N2+ + O',
   fr={'Arrhenius', A=7.2e13,  n=0.00, C=35500.0},
}

-- Electron-impact ionization reactions
Reaction{'O + e- <=> O+ + e- + e-',
   fr={'Park', A=3.9e33, n=-3.78, C=158500.0, s=0.0},
}

Reaction{'N + e- <=> N+ + e- + e-',
   fr={'Park', A=2.5e34, n=-3.82, C=168600.0, s=0.0},
}

-- Radiative recombination reactions
-- Dan had these commented out for some reason. For the sake of comparison we'll leave them out too.
-- Reaction{'O+ + e- <=> O',
--     fr={'Park', A=1.07e11, n=-0.52, C=0.0, s=0.0},
-- }

-- Reaction{'N+ + e- <=> N',
--     fr={'Park', A=1.52e11, n=-0.48, C=0.0, s=0.0},
-- }
