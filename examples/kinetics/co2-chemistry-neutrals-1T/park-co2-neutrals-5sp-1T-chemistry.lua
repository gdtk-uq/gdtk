--
-- This chemical kinetic system provides
-- a HACKED,PARTIAL reaction scheme for the Martian atmosphere.
--
-- Reference: Park, Howe, Jaffe and Candler - 
-- Review of Chemical-Kinetic Problems of Future NASA Missions, II: Mars Entries 
-- https://arc.aiaa.org/doi/pdf/10.2514/3.496
--
-- Author of the original unhacked file: Jens Kunze
-- Date: 2020-Jan-16
-- Place: Brisbane, Australia
-- This reduced version: Rowan G. and Peter J.
--
-- Species are CO2, O2, CO, O, C 
--
-- Reaction rates in table 2 in the reference paper are
-- given in the format 3.7^14. In the text, however,
-- they are used as 3.7e+14!
--
Config{
   tempLimits = {lower=300.0, upper=15000.0}
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
    'CO + CO2 <=> C + O + CO2',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r36'
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
    'CO2 + O2 <=> CO + O + O2',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r54'
}

Reaction{
    'CO2 + CO <=> CO + O + CO',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r56'
}

Reaction{
    'CO2 + CO2 <=> CO + O + CO2',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r58'
}

Reaction{
    'CO + O <=> C + O2',
    fr={'Arrhenius', A=3.9e+13, n=-0.18, C=69200},
    label='r62'
}

Reaction{
    'CO2 + O <=> O2 + CO',
    fr={'Arrhenius', A=2.1e+13, n=0.00, C=27800},
    label='r68'
}
