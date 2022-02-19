-- Author: Rowan J. Gollan
-- Date: 2021-06-12
--
-- Reference:
--   Ebrahim and Hornung (1973)
--   Nonequilibrium Nozzle Expansions of Carbon Dioxide from a
--   High-Enthalpy Reservoir
--   AIAA Journal, 11(10), pp. 1369--1370
--
-- From Table 1.

R_univ = 1.987e-3 -- kcal/(mol.K) [from Wikipedia]

-- Dissociation --

Reaction{
   'CO2 + M <=> CO + O + M',
   fr={'Park', A=2.88e11,  n=0.5,  C=74.8/R_univ, s=0.5},
   label='r1'
}

Reaction{
   'CO + M <=> C + O + M',
   fr={'Park', A=8.79e29, n=-3.52, C=255.76/R_univ, s=0.5},
   efficiencies={CO=0, O=0}, -- because these are excluded as M, they have separate rates.
   label='r2'
}

Reaction{
   'CO + CO <=> C + O + CO',
   fr={'Park', A=1.76e30, n=-3.52, C=255.76/R_univ, s=0.5},
   label='r3'
}

Reaction{
   'CO + O <=> C + O + O',
   fr={'Park', A=1.29e31, n=-3.52, C=255.76/R_univ, s=0.5},
   label='r4'
}

Reaction{
   'O2 + M <=> 2O + M',
   fr={'Park', A=2.55e18, n=-1.0, C=118.7/R_univ, s=0.5},
   efficiencies={O=0, O2=0}, -- because these are excluded as M, they have separate rates.
   label='r5'
}

Reaction{
   'O2 + O2 <=> 2O + O2',
   fr={'Park', A=2.75e19, n=-1.0, C=118.7/R_univ, s=0.5},
   label='r6'
}

Reaction{
   'O2 + O <=> 3O',
   fr={'Park', A=2.1e18, n=-1.0, C=117.96/R_univ, s=0.5},
   label='r7'
}

-- Carbon exchange --

Reaction{
   'CO + CO <=> CO2 + C',
   fr={'Arrhenius', A=2.33e9, n=0.5, C=130.5/R_univ},
   label='r8'
}

Reaction{
   'O + CO <=> C + O2',
   fr={'Arrhenius', A=2.73e11, n=0.5, C=138.1/R_univ},
   label='r9'
}

Reaction{
   'CO + O2 <=> CO2 + O',
   fr={'Arrhenius', A=1.6e13, n=0, C=41.0/R_univ},
   label='r10'
}
