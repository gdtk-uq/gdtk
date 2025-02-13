--
-- This chemical kinetic system provides
-- a full reaction scheme for the Martian atmosphere.
--
-- Reference: Park, Howe, Jaffe and Candler - 
-- Review of Chemical-Kinetic Problems of Future NASA Missions, II: Mars Entries 
-- https://arc.aiaa.org/doi/pdf/10.2514/3.496
--
-- Author: Jens Kunze
-- Date: 2020-Jan-16
-- Place: BRisbane, Australia
--
-- Species are CO2, N2, Ar, O2, CO, O, C, N, CN, C2, NO, NCO, 
--              NO+, O2+, CO+, O+, C+, e-
--
-- Ar+, N+, C2+, N2+ and CN+ are possible but their concentrations
-- are insignificantly small according to Park et al
--
-- Reaction rates in table 2 in the reference paper are
-- given in the format 3.7^14. In the text, however,
-- they are used as 3.7e+14!
--
-- Comment NO_SPECIES=12 if you want to use the full reaction scheme
-- with ions
----------------------------------------------------------------

-- RJG, 2025-02-13
-- Mechanism reduced to neutrals only

Config{
   tempLimits = {lower=300.0, upper=15000.0}
}


Reaction{
    'CO2 + C <=> CO + O + C',
    fr={'Arrhenius', A=1.4e+22, n=-1.50, C=63275},
    label='r49'
}

Reaction{
    'CO2 + O <=> CO + O + O',
    fr={'Arrhenius', A=1.4e+22, n=-1.50, C=63275},
    label='r51'
}


Reaction{
    'CO + C <=> C + O + C',
    fr={'Arrhenius', A=3.4e+20, n=-1.00, C=129000},
    label='r27'
}

Reaction{
    'CO + O <=> C + O + O',
    fr={'Arrhenius', A=3.4e+20, n=-1.00, C=129000},
    label='r29'
}

Reaction{
    'CO + O2 <=> C + O + O2',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r32'
}

Reaction{
    'CO + CO <=> C + O + CO',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r34'
}


Reaction{
    'CO + CO <=> C + O + CO',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r34'
}

Reaction{
    'O2 + C <=> O + O + C',
    fr={'Arrhenius', A=1.0e+22, n=-1.50, C=59750},
    label='r15'
}

Reaction{
    'O2 + O <=> O + O + O',
    fr={'Arrhenius', A=1.0e+22, n=-1.50, C=59750},
    label='r17'
}

Reaction{
    'O2 + O2 <=> O + O + O2',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r20'
}

Reaction{
    'O2 + CO <=> O + O + CO',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r22'
}

Reaction{
    'O2 + CO2 <=> O + O + CO2',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r24'
}

Reaction{
    'CO2 + O <=> O2 + CO',
    fr={'Arrhenius', A=2.1e+13, n=0.00, C=27800},
    label='r68'
}
