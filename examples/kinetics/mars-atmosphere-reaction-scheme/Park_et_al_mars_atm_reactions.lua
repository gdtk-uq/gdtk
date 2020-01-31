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

NO_SPECIES = 12

-- Currently two schemes are available in this file using 
-- the WITH_IONIZATION flag. If 'false' then 72
-- reactions are used. If 'true' the full 84 reactions
-- are considered. If other schemes are required this can
-- be configured in the same manner using the 
-- 'select_reactions_by_label' function as demonstrated at
-- the end of the file.

WITH_IONIZATION = true

if NO_SPECIES == 12 then
   WITH_IONIZATION = false
end

-- We know that some of the species in the efficiencies list do NOT
-- apply when we don't consider ionization.  We choose to suppress
-- those warnings.

-- SUPPRESS_WARNINGS = true

-- Temperature limits.
-- These temperature limits are used to limit the evaluation of
-- the rate coefficients. Park et al state that up to 15000K
-- the radiation is evaluated adequately, so it may not be wise
-- to evaluate these rate coefficients for higher temperatures. Also, the
-- rate values do get meaningless at very low temperatures.
-- For this reason, we suggest the limits we set here.
-- WARNING: We have experienced a flow solver/chemistry coupling
-- issue for strong expansions when the lower limit is set too low.
-- 

Config{
   tempLimits = {lower=300.0, upper=15000.0}
}

-- Dissociation Reactions
Reaction{
    'C2 + M <=> C + C + M',
    fr={'Arrhenius', A=3.7e+14, n=0.00, C=69900},
    label='r1'
}

Reaction{
    'N2 + Ar <=> N + N + Ar',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r2'
}

Reaction{
    'N2 + C <=> N + N + C',
    fr={'Arrhenius', A=3.0e+22, n=-1.60, C=113200},
    label='r3'
}

Reaction{
    'N2 + N <=> N + N + N',
    fr={'Arrhenius', A=3.0e+22, n=-1.60, C=113200},
    label='r4'
}

Reaction{
    'N2 + O <=> N + N + O',
    fr={'Arrhenius', A=3.0e+22, n=-1.60, C=113200},
    label='r5'
}

Reaction{
    'N2 + C2 <=> N + N + C2',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r6'
}

Reaction{
    'N2 + N2 <=> N + N + N2',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r7'
}

Reaction{
    'N2 + O2 <=> N + N + O2',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r8'
}

Reaction{
    'N2 + CN <=> N + N + CN',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r9'
}

Reaction{
    'N2 + CO <=> N + N + CO',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r10'
}

Reaction{
    'N2 + NO <=> N + N + NO',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r11'
}

Reaction{
    'N2 + CO2 <=> N + N + CO2',
    fr={'Arrhenius', A=7.0e+21, n=-1.60, C=113200},
    label='r12'
}

Reaction{
    'N2 + e- <=> N + N + e-',
    fr={'Arrhenius', A=1.25e+25, n=-1.60, C=113200},
    label='r13'
}

Reaction{
    'O2 + Ar <=> O + O + Ar',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r14'
}

Reaction{
    'O2 + C <=> O + O + C',
    fr={'Arrhenius', A=1.0e+22, n=-1.50, C=59750},
    label='r15'
}

Reaction{
    'O2 + N <=> O + O + N',
    fr={'Arrhenius', A=1.0e+22, n=-1.50, C=59750},
    label='r16'
}

Reaction{
    'O2 + O <=> O + O + O',
    fr={'Arrhenius', A=1.0e+22, n=-1.50, C=59750},
    label='r17'
}

Reaction{
    'O2 + C2 <=> O + O + C2',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r18'
}

Reaction{
    'O2 + N2 <=> O + O + N2',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r19'
}

Reaction{
    'O2 + O2 <=> O + O + O2',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r20'
}

Reaction{
    'O2 + CN <=> O + O + CN',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r21'
}

Reaction{
    'O2 + CO <=> O + O + CO',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r22'
}

Reaction{
    'O2 + NO <=> O + O + NO',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r23'
}

Reaction{
    'O2 + CO2 <=> O + O + CO2',
    fr={'Arrhenius', A=2.0e+21, n=-1.50, C=59750},
    label='r24'
}

Reaction{
    'CN + M <=> C + N + M',
    fr={'Arrhenius', A=2.5e+14, n=0.00, C=71000},
    label='r25',
}

Reaction{
    'CO + Ar <=> C + O + Ar',
    fr={'Arrhenius', A=2.3e+19, n=-1.00, C=129000},
    label='r26'
}

Reaction{
    'CO + C <=> C + O + C',
    fr={'Arrhenius', A=3.4e+20, n=-1.00, C=129000},
    label='r27'
}

Reaction{
    'CO + N <=> C + O + N',
    fr={'Arrhenius', A=3.4e+20, n=-1.00, C=129000},
    label='r28'
}

Reaction{
    'CO + O <=> C + O + O',
    fr={'Arrhenius', A=3.4e+20, n=-1.00, C=129000},
    label='r29'
}

Reaction{
    'CO + C2 <=> C + O + C2',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r30'
}

Reaction{
    'CO + N2 <=> C + O + N2',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r31'
}

Reaction{
    'CO + O2 <=> C + O + O2',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r32'
}

Reaction{
    'CO + CN <=> C + O + CN',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r33'
}

Reaction{
    'CO + CO <=> C + O + CO',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r34'
}

Reaction{
    'CO + NO <=> C + O + NO',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r35'
}

Reaction{
    'CO + CO2 <=> C + O + CO2',
    fr={'Arrhenius', A=2.3e+20, n=-1.00, C=129000},
    label='r36'
}

Reaction{
    'NO + Ar <=> N + O + Ar',
    fr={'Arrhenius', A=5.0e+15, n=0.00, C=75500},
    label='r37'
}

Reaction{
    'NO + C <=> N + O + C',
    fr={'Arrhenius', A=1.1e+17, n=0.00, C=75500},
    label='r38'
}

Reaction{
    'NO + N <=> N + O + N',
    fr={'Arrhenius', A=1.1e+17, n=0.00, C=75500},
    label='r39'
}

Reaction{
    'NO + O <=> N + O + O',
    fr={'Arrhenius', A=1.1e+17, n=0.00, C=75500},
    label='r40'
}

Reaction{
    'NO + C2 <=> N + O + C2',
    fr={'Arrhenius', A=5.0e+15, n=0.00, C=75500},
    label='r41'
}

Reaction{
    'NO + N2 <=> N + O + N2',
    fr={'Arrhenius', A=5.0e+15, n=0.00, C=75500},
    label='r42'
}

Reaction{
    'NO + O2 <=> N + O + O2',
    fr={'Arrhenius', A=5.0e+15, n=0.00, C=75500},
    label='r43'
}

Reaction{
    'NO + CN <=> N + O + CN',
    fr={'Arrhenius', A=5.0e+15, n=0.00, C=75500},
    label='r44'
}

Reaction{
    'NO + CO <=> N + O + CO',
    fr={'Arrhenius', A=5.0e+15, n=0.00, C=75500},
    label='r45'
}

Reaction{
    'NO + NO <=> N + O + NO',
    fr={'Arrhenius', A=1.1e+17, n=0.00, C=75500},
    label='r46'
}

Reaction{
    'NO + CO2 <=> N + O + CO2',
    fr={'Arrhenius', A=1.1e+17, n=0.00, C=75500},
    label='r47'
}

Reaction{
    'CO2 + Ar <=> CO + O + Ar',
    fr={'Arrhenius', A=6.9e+20, n=-1.50, C=63275},
    label='r48'
}

Reaction{
    'CO2 + C <=> CO + O + C',
    fr={'Arrhenius', A=1.4e+22, n=-1.50, C=63275},
    label='r49'
}

Reaction{
    'CO2 + N <=> CO + O + N',
    fr={'Arrhenius', A=1.4e+22, n=-1.50, C=63275},
    label='r50'
}

Reaction{
    'CO2 + O <=> CO + O + O',
    fr={'Arrhenius', A=1.4e+22, n=-1.50, C=63275},
    label='r51'
}

Reaction{
    'CO2 + C2 <=> CO + O + C2',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r52'
}

Reaction{
    'CO2 + N2 <=> CO + O + N2',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r53'
}

Reaction{
    'CO2 + O2 <=> CO + O + O2',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r54'
}

Reaction{
    'CO2 + CN <=> CO + O + CN',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r55'
}

Reaction{
    'CO2 + CO <=> CO + O + CO',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r56'
}

Reaction{
    'CO2 + NO <=> CO + O + NO',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r57'
}

Reaction{
    'CO2 + CO2 <=> CO + O + CO2',
    fr={'Arrhenius', A=6.9e+21, n=-1.50, C=63275},
    label='r58'
}

Reaction{
    'NCO + M <=> CO + N + M',
    fr={'Arrhenius', A=6.3e+16, n=-0.50, C=24000},
    label='r59'
}

-- Neutral exchange reactions
Reaction{
    'NO + O <=> N + O2',
    fr={'Arrhenius', A=8.4e+12, n=0.00, C=19450},
    label='r60'
}

Reaction{
    'N2 + O <=> NO + N',
    fr={'Arrhenius', A=6.4e+17, n=-1.00, C=38370},
    label='r61'
}

Reaction{
    'CO + O <=> C + O2',
    fr={'Arrhenius', A=3.9e+13, n=-0.18, C=69200},
    label='r62'
}

Reaction{
    'CO + C <=> C2 + O',
    fr={'Arrhenius', A=2.0e+17, n=-1.00, C=58000},
    label='r63'
}

Reaction{
    'CO + N <=> CN + O',
    fr={'Arrhenius', A=1.0e+14, n=0.00, C=38600},
    label='r64'
}

Reaction{
    'N2 + C <=> CN + N',
    fr={'Arrhenius', A=1.1e+14, n=-0.11, C=23200},
    label='r65'
}

Reaction{
    'CN + O <=> NO + C',
    fr={'Arrhenius', A=1.6e+13, n=0.10, C=14600},
    label='r66'
}

Reaction{
    'CN + C <=> C2 + N',
    fr={'Arrhenius', A=5.0e+13, n=0.00, C=13000},
    label='r67'
}

Reaction{
    'CO2 + O <=> O2 + CO',
    fr={'Arrhenius', A=2.1e+13, n=0.00, C=27800},
    label='r68'
}

Reaction{
    'CN + O2 <=> NCO + O',
    fr={'Arrhenius', A=6.6e+12, n=0.00, C=-200},
    label='r69'
}

Reaction{
    'CN + CO2 <=> NCO + CO',
    fr={'Arrhenius', A=4.0e+14, n=0.00, C=19200},
    label='r70'
}

Reaction{
    'CN + NO <=> NCO + N',
    fr={'Arrhenius', A=1.0e+14, n=0.00, C=21200},
    label='r71'
}

Reaction{
    'CO + NO <=> NCO + O',
    fr={'Arrhenius', A=3.8e+17, n=-0.873, C=51600},
    label='r72'
}

Reaction{
    'CN + CO <=> NCO + C',
    fr={'Arrhenius', A=1.5e+16, n=-0.487, C=65800},
    label='r73'
}

-- Associative ionization reactions
Reaction{
    'N + O <=> NO+ + e-',
    fr={'Arrhenius', A=8.8e+8, n=1.00, C=31900},
    label='r74'
}

Reaction{
    'O + O <=> O2+ + e-',
    fr={'Arrhenius', A=7.1e+2, n=2.70, C=80600},
    label='r75'
}

Reaction{
    'C + O <=> CO+ + e-',
    fr={'Arrhenius', A=8.8e+8, n=1.00, C=33100},
    label='r76'
}

-- Charge exchange reactions
Reaction{
    'NO+ + C <=> NO + C+',
    fr={'Arrhenius', A=1.0e+13, n=0.00, C=23200},
    label='r77'
}

Reaction{
    'O2+ + O <=> O+ + O2',
    fr={'Arrhenius', A=4.0e+12, n=-0.09, C=18000},
    label='r78'
}

Reaction{
    'NO+ + N <=> O+ + N2',
    fr={'Arrhenius', A=3.4e+13, n=-1.08, C=12800},
    label='r79'
}

Reaction{
    'NO+ + O <=> O2+ + N',
    fr={'Arrhenius', A=7.2e+12, n=0.29, C=48600},
    label='r80'
}

Reaction{
    'CO + C+ <=> CO+ + C',
    fr={'Arrhenius', A=1.0e+13, n=0.00, C=31400},
    label='r81'
}

Reaction{
    'O2 + C+ <=> O2+ + C',
    fr={'Arrhenius', A=1.0e+13, n=0.00, C=9400},
    label='r82'
}

-- Electron-impact ionization reactions
Reaction{
    'C + e- <=> C+ + e- + e-',
    fr={'Arrhenius', A=3.9e+33, n=-3.78, C=130700},
    label='r83'
}

Reaction{
    'O + e- <=> O+ + e- + e-',
    fr={'Arrhenius', A=3.9e+33, n=-3.78, C=158500},
    label='r84'
}

--Radiative recombination reactions
--commented for the moment because there is no way to treat the hv yet
--
--Reaction{
--    'O+ + e- <=> O + hv',
--    fr={'Arrhenius', A=1.07e+11, n=-0.52, C=0},
--    label='r85'
--}
--
--Reaction{
--    'C+ + e- <=> C + hv',
--    fr={'Arrhenius', A=2.02e+11, n=-0.46, C=0},
--    label='r86'
--}

if NO_SPECIES == 12 then
    list_of_reactions = {}
    for i=1, 73 do
        if i ~= 13 then
            list_of_reactions[#list_of_reactions+1] = 'r' .. i
        end
    end
   selectOnlyReactionsWithLabel(list_of_reactions)
end
