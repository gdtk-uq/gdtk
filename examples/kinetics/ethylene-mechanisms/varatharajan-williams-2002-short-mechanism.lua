-- Author: Rowan J. Gollan
-- Date: 2023-06-11
--
-- This file represents Table 1: Short Mechanism for C2H4 detonation
-- in Varatharajan and Williams (2002).
--
-- Reference:
-- Varatharajan and Williams (2002)
-- Ethylene Ignition and Detonation Chemistry, Part 2:
-- Ignition Histories and Reduced Mechanims
-- Journal of Propulsion and Power, 18(2), pp. 352--362
--
-- Notes:
-- 1. Table 1 in V-W (2002) gives activation energy, E,
--    in kJ/mol. So we'll convert to an activation temperature, C,
--    as C = E/R
--
-- 2. All reactions are listed in a forward direction only,
--    but some are couplets indicating forward and reverse.
--    For example:
--    1: H + O2 --> OH + O
--    2: OH + O --> H + O2
--    It is tempting to list these as one reaction with a
--    a backwards rate constant, however, the units for
--    A may not be correct depending on the colliders.
--
--    Instead, we'll list all reactions as per Table 1.
--    The backwards rate constant will be set to zero
--    to indicate that the reaction is not reversible.

R = 8.314e-3 -- kJ/(K.mol)

irreversible = {'Arrhenius', A=0.0, n=0.0, C=0.0}

Config{
   odeStep = {method='alpha-qss'}
}

-----
-- Hydrogen-oxygen chain
-----

Reaction{
   'H + O2 => OH + O',
   fr={'Arrhenius', A=3.52e16, n=-0.70, C=71.4/R},
   br=irreversible,
   label='r1'
}

Reaction{
   'OH + O => H + O2',
   fr={'Arrhenius', A=1.15e14, n=-0.32, C=-0.7/R},
   br=irreversible,
   label='r2'
}

Reaction{
   'OH + H2 => H2O + H',
   fr={'Arrhenius', A=1.17e9, n=1.30, C=15.2/R},
   br=irreversible,
   label='r3'
}

Reaction{
   'O + H2O => 2OH',
   fr={'Arrhenius', A=7.60e0, n=3.84, C=53.5/R},
   br=irreversible,
   label='r4'
}

Reaction{
   '2OH => O + H2O',
   fr={'Arrhenius', A=2.45e-1, n=3.97, C=-19.0/R},
   br=irreversible,
   label='r5'
}

-----
-- Hydroperoxyl formation and consumption
-----

Reaction{
   'H + O2 + M => HO2 + M',
   fr={'Arrhenius', A=2.60e19, n=-1.20, C=0.0/R},
   br=irreversible,
   efficiencies={ O2=0.3, H2O=7.0, CO=0.75, CO2=1.5, C2H6=1.5, Ar=0.5 },
   label='r6'
}

Reaction{
   'HO2 + H => 2OH',
   fr={'Arrhenius', A=1.70e14, n=0.00, C=3.7/R},
   br=irreversible,
   label='r7'
}

Reaction{
   'HO2 + H => H2 + O2',
   fr={'Arrhenius', A=4.28e13, n=0.00, C=5.9/R},
   br=irreversible,
   label='r8'
}

Reaction{
   'HO2 + OH => H2O + O2',
   fr={'Arrhenius', A=2.89e13, n=0.00, C=-2.1/R},
   br=irreversible,
   label='r9'
}

-----
-- Hydrogen peroxide formation and consumption
-----

Reaction{
   '2HO2 => H2O2 + O2',
   fr={'Arrhenius', A=3.02e12, n=0.00, C=5.8/R},
   br=irreversible,
   label='r10'
}

Reaction{
   'H2O2 (+ M) => 2OH (+ M)',
   fr={'pressure dependent',
       k0=   {A=7.94e24, n=-2.21, C=212.0/R},
       kInf= {A=2.55e20, n=-1.68, C=219.1/R},
       Troe= {a=0.735, T1=1756, T2=5182, T3=94}
   },
   br=irreversible,
   efficiencies={ H2=2.0, H2O=6.0, CO=1.5, CO2=2.0, Ar=0.7 },
   label='r11'
}

-----
-- Direct recombination
-----

Reaction{
   'H + OH + M => H2O + M',
   fr={'Arrhenius', A=2.20e22, n=-2.00, C=0.0/R},
   br=irreversible,
   efficiencies={ O2=0.3, H2O=7.0, CO=0.75, CO2=1.55, C2H6=1.5, Ar=0.5 },
   label='r12'
}

Reaction{
   'H2O + M => H + OH + M',
   fr={'Arrhenius', A=2.18e23, n=-1.93, C=499.0/R},
   br=irreversible,
   efficiencies={ O2=0.3, H2O=7.0, CO=0.75, CO2=1.55, C2H6=1.5, Ar=0.5 },
   label='r13'
}

-----
-- Carbon monoxide reactions
-----

Reaction{
   'CO + OH => CO2 + H',
   fr={'Arrhenius', A=4.40e6, n=1.50, C=-3.1/R},
   br=irreversible,
   label='r14'
}

Reaction{
   'CO2 + H => CO + OH',
   fr={'Arrhenius', A=4.97e8, n=1.50, C=89.7/R},
   br=irreversible,
   label='r15'
}

-----
-- Initiation and fuel consumption
-----

Reaction{
   'C2H4 + O2 => C2H3 + HO2',
   fr={'Arrhenius', A=4.22e13, n=0.0, C=241.0/R},
   br=irreversible,
   label='r16'
}

Reaction{
   'C2H4 + OH => C2H3 + H2O',
   fr={'Arrhenius', A=5.53e05, n=2.31, C=12.4/R},
   br=irreversible,
   label='r17'
}

Reaction{
   'C2H4 + O => CH3 + CHO',
   fr={'Arrhenius', A=2.25e06, n=2.08, C=0.0/R},
   br=irreversible,
   label='r18'
}

Reaction{
   'C2H4 + O => CH2CHO + H',
   fr={'Arrhenius', A=1.21e06, n=2.08, C=0.0/R},
   br=irreversible,
   label='r19'
}

Reaction{
   'C2H4 + HO2 => C2H4O + OH',
   fr={'Arrhenius', A=2.23e12, n=0.00, C=71.9/R},
   br=irreversible,
   label='r20'
}

Reaction{
   'C2H4 + H => C2H3 + H2',
   fr={'Arrhenius', A=4.49e07, n=2.12, C=55.9/R},
   br=irreversible,
   label='r21'
}

Reaction{
   'C2H4 + H (+ M) => C2H5 (+ M)',
   fr={'pressure dependent',
       k0=   {A=1.90e35, n=-5.57, C=21.1/R},
       kInf= {A=1.08e12, n=0.45, C=7.6/R},
       Troe= {a=(1 - 0.832), T3=1203}
   },
   br=irreversible,
   label='r22'
}

Reaction{
   'C2H4 + M => C2H3 + H + M',
   fr={'Arrhenius', A=2.60e17, n=0.00, C=404.0/R},
   br=irreversible,
   efficiencies={ H2=2.0, H2O=6.0, CO=1.5, CO2=2.0, Ar=0.7 },
   label='r23'
}

-----
-- Vinyl, methyl, vinoxy, and ethyl consumption
-----

Reaction{
   'C2H3 + H => C2H2 + H2',
   fr={'Arrhenius', A=1.21e13, n=0.00, C=0.0/R},
   br=irreversible,
   label='r24'
}

Reaction{
   'C2H3 + O2 => CH2O + CHO',
   fr={'Arrhenius', A=1.70e29, n=-5.31, C=27.2/R},
   br=irreversible,
   label='r25'
}

Reaction{
   'C2H3 + O2 => CH2CHO + O',
   fr={'Arrhenius', A=7.00e14, n=-0.61, C=22.0/R},
   br=irreversible,
   label='r26'
}

Reaction{
   'CH3 + O2 => CH2O + OH',
   fr={'Arrhenius', A=3.30e11, n=0.00, C=37.4/R},
   br=irreversible,
   label='r27'
}

Reaction{
   'CH3 + O => CH2O + H',
   fr={'Arrhenius', A=8.43e13, n=0.00, C=0.0/R},
   br=irreversible,
   label='r28'
}

Reaction{
   'CH2CHO =>  CH2CO + H',
   fr={'Arrhenius', A=1.05e37, n=-7.19, C=186.0/R},
   br=irreversible,
   label='r29'
}

Reaction{
   'C2H5 + O2 => C2H4 + HO2',
   fr={'Arrhenius', A=2.00e12, n=0.00, C=20.9/R},
   br=irreversible,
   label='r30'
}

Reaction{
   'C2H5 (+ M) => C2H4 + H (+ M)',
   fr={'pressure dependent',
       k0=   {A=3.99e33, n=-4.99, C=167.4/R},
       kInf= {A=1.11e10, n=1.04, C=153.8/R},
       Troe= {a=(1 - 0.832), T3=1203}
   },
   br=irreversible,
   label='r31'
}

-----
-- Ketene, formaldehyde, formyl, acetylene, and ethylene oxide consumption
-----

Reaction{
   'CH2CO + H => CH3 + CO',
   fr={'Arrhenius', A=1.11e07, n=2.00, C=8.4/R},
   br=irreversible,
   label='r32'
}

Reaction{
   'CH2O + OH => CHO + H2O',
   fr={'Arrhenius', A=3.90e10, n=0.89, C=1.7/R},
   br=irreversible,
   label='r33'
}

Reaction{
   'CHO + M => CO + H + M',
   fr={'Arrhenius', A=1.86e17, n=-1.00, C=71.1/R},
   br=irreversible,
   efficiencies={ CO=2.5, CO2=2.5, H2=1.9, H2O=12.0 },
   label='r34'
}

Reaction{
   'CHO + O2 => CO + HO2',
   fr={'Arrhenius', A=3.00e12, n=0.00, C=0.0/R},
   br=irreversible,
   label='r35'
}

Reaction{
   'CHO + H => CO + H2',
   fr={'Arrhenius', A=1.00e14, n=0.00, C=0.0/R},
   br=irreversible,
   label='r36'
}

Reaction{
   'C2H2 + OH => CH2CO + H',
   fr={'Arrhenius', A=1.90e07, n=1.70, C=4.2/R},
   br=irreversible,
   label='r37'
}

Reaction{
   'C2H4O + HO2 => CH3 + CO + H2O2',
   fr={'Arrhenius', A=4.00e12, n=0.00, C=71.2/R},
   br=irreversible,
   label='r38'
}




