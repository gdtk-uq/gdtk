-- Author: Rowan J. Gollan
-- Date: 29-Mar-2009
-- Place: Poquoson, Virginia, USA
--
-- Adapted from Python file: rogers_schexnayder.py
--
-- This file provides a reaction scheme for
-- hydrogen combustion in air.
-- NOTE: This scheme does not include carbonaceous compounds
-- or Argon (or any of the associated reactions).
--
-- Reference:
-- Rogers, R.C. and Schexnayder, Jr., C.J. (1981)
-- Chemical Kinetic Analysis of Hydroden-Air
-- Ignition and Reaction Times
-- NASA Technical Paper 1856
--
-- Species used: O, O2, N, N2, H, H2, H2O, HO2, OH, NO, NO2, HNO2, HNO3, O3, H2O2, HNO
--
-- Updated 2016-07-31
--    Updated for eilmer4

options = {}
options.H2_O2_only = false

Config{
   odeStep = {method='alpha-qss'},
   tightTempCoupling = true
}

Reaction{
   'O2 + M <=> O + O + M',
   fr={"Arrhenius", A=0.72e19, n=-1.0, C=59340.0},
   efficiencies={O2=4.0, O=10.0, H2O=2.0},
   label='r1'
}

Reaction{
   'M + H2 <=> H + H + M',
   fr={'Arrhenius', A=0.55e19, n=-1.0, C=51987.0},
   efficiencies={H=5.0, H2=2.0, H2O=8.0},
   label='r2'
}

Reaction{
   'M + H2O <=> H + OH + M',
   fr={'Arrhenius', A=0.52e22, n=-1.5, C=59386.0},
   efficiencies={H2O=6.0},
   label='r3'
}

Reaction{
   'H + O2 + M <=> HO2 + M',
   fr={'Arrhenius', A=0.23e16, n=0.0, C=-403.0},
   efficiencies={H2=2.0, H2O=13.0},
   label='r4'
}

Reaction{
   'M + NO2 <=> NO + O + M',
   fr={'Arrhenius', A=0.11e17, n=0.0, C=32710.0},
   label='r5'
}

Reaction{
   'M + NO <=> N + O + M',
   fr={'Arrhenius', A=0.41e19, n=-1.0, C=75330.0},
   label='r6'
}

-- not included, ignoring the carbonaceous compounds
r7 = 'M + O + CO <=> CO2 + M'

Reaction{
   'M + H + NO <=> HNO + M',
   fr={'Arrhenius', A=0.54e16, n=0.0, C=-300.0},
   efficiencies={H2O=3.0},
   label='r8'
}

Reaction{
   'M + H2O2 <=> OH + OH + M',
   fr={'Arrhenius', A=0.12e18, n=0.0, C=22899.0},
   efficiencies={H2O=6.0},
   label='r9'
}

Reaction{
   'M + OH + NO <=> HNO2 + M',
   fr={'Arrhenius', A=0.80e16, n=0.0, C=-1000.0},
   label='r10'
}

Reaction{
   'M + OH + NO2 <=> HNO3 + M',
   fr={'Arrhenius', A=0.13e17, n=0.0, C=-1107.0},
   label='r11'
}

Reaction{
   'M + O3 <=> O2 + O + M',
   fr={'Arrhenius', A=0.13e22, n=-2.0, C=12800.0},
   efficiencies={O2=1.5},
   label='r12'
}

r13 = 'M + HCO <=> CO + H + M'

Reaction{
   'M + O + H <=> OH + M',
   fr={'Arrhenius', A=0.71e19, n=-1.0, C=0.0},
   label='r14'
}

Reaction{
   'H2O + O <=> OH + OH',
   fr={'Arrhenius', A=0.58e14, n=0.0, C=9059.0},
   label='r15'
}

Reaction{
   'H2 + OH <=> H2O + H',
   fr={'Arrhenius', A=0.20e14, n=0.0, C=2600.0},
   label='r16'
}

Reaction{
   'O2 + H <=> OH + O',
   fr={'Arrhenius', A=0.22e15, n=0.0, C=8455.0},
   label='r17'
}

Reaction{
   'H2 + O <=> OH + H',
   fr={'Arrhenius', A=0.75e14, n=0.0, C=5586.0},
   label='r18'
}

Reaction{
   'H2 + O2 <=> OH + OH',
   fr={'Arrhenius', A=0.10e14, n=0.0, C=21641.0},
   label='r19'
}

Reaction{
   'H + HO2 <=> H2 + O2',
   fr={'Arrhenius', A=0.24e14, n=0.0, C=350.0},
   label='r20'
}

Reaction{
   'H2 + O2 <=> H2O + O',
   fr={'Arrhenius', A=0.41e14, n=0.0, C=25400.0},
   label='r21'
}

Reaction{
   'H + HO2 <=> OH + OH',
   fr={'Arrhenius', A=0.24e15, n=0.0, C=950.0},
   label='r22'
}

Reaction{
   'H2O + O <=> H + HO2',
   fr={'Arrhenius', A=0.58e12, n=0.5, C=28686.0},
   label='r23'
}

Reaction{
   'O + HO2 <=> OH + O2',
   fr={'Arrhenius', A=0.5e14, n=0.0, C=503.0},
   label='r24'
}

Reaction{
   'OH + HO2 <=> O2 + H2O',
   fr={'Arrhenius', A=0.3e14, n=0.0, C=0.0},
   label='r25'
}

Reaction{
   'H2 + HO2 <=> H2O + OH',
   fr={'Arrhenius', A=0.2e14, n=0.0, C=12582.0},
   label='r26'
}

Reaction{
   'HO2 + H2 <=> H + H2O2',
   fr={'Arrhenius', A=0.73e12, n=0.0, C=9400.0},
   label='r27'
}

Reaction{
   'H2O2 + H <=> OH + H2O',
   fr={'Arrhenius', A=0.32e15, n=0.0, C=4504.0},
   label='r28'
}

Reaction{
   'HO2 + OH <=> O + H2O2',
   fr={'Arrhenius', A=0.52e11, n=0.5, C=10600.0},
   label='r29'
}

Reaction{
   'HO2 + H2O <=> OH + H2O2',
   fr={'Arrhenius', A=0.28e14, n=0.0, C=16500.0},
   label='r30'
}

Reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=0.2e13, n=0.0, C=0.0},
   label='r31'
}

Reaction{
   'O + O3 <=> O2 + O2',
   fr={'Arrhenius', A=0.10e14, n=0.0, C=2411.0},
   label='r32'
}

Reaction{
   'O3 + NO <=> NO2 + O2',
   fr={'Arrhenius', A=0.54e12, n=0.0, C=1200.0},
   label='r33'
}

Reaction{
   'O3 + H <=> OH + O2',
   fr={'Arrhenius', A=0.70e14, n=0.0, C=560.0},
   label='r34'
}

Reaction{
   'O3 + OH <=> O2 + HO2',
   fr={'Arrhenius', A=0.90e12, n=0.0, C=1000.0},
   label='r35'
}

Reaction{
   'O + N2 <=> NO + N',
   fr={'Arrhenius', A=0.50e14, n=0.0, C=37940.0},
   label='r36'
}

Reaction{
   'H + NO <=> OH + N',
   fr={'Arrhenius', A=0.17e15, n=0.0, C=24500.0},
   label='r37'
}

Reaction{
   'O + NO <=> O2 + N',
   fr={'Arrhenius', A=0.15e10, n=1.0, C=19500.0},
   label='r38'
}

Reaction{
   'NO2 + H <=> NO + OH',
   fr={'Arrhenius', A=0.35e15, n=0.0, C=740.0},
   label='r39'
}

Reaction{
   'NO2 + O <=> NO + O2',
   fr={'Arrhenius', A=0.10e14, n=0.0, C=302.0},
   label='r40'
}

Reaction{
   'NO2 + H2 <=> HNO2 + H',
   fr={'Arrhenius', A=0.24e14, n=0.0, C=14595.0},
   label='r41'
}

Reaction{
   'HO2 + NO <=> NO2 + OH',
   fr={'Arrhenius', A=0.30e13, n=0.5, C=1208.0},
   label='r42'
}

Reaction{
   'NO2 + H2O <=> HNO2 + OH',
   fr={'Arrhenius', A=0.32e13, n=0.0, C=22000.0},
   label='r43'
}

Reaction{
   'NO2 + OH <=> HNO2 + O',
   fr={'Arrhenius', A=0.21e13, n=0.0, C=12580.0},
   label='r44'
}

r45 = 'CO + OH <=> CO2 + H'
r46 = 'CO2 + O <=> O2 + CO'
r47 = 'H2O + CO <=> HCO + OH'
r48 = 'OH + CO <=> HCO + O'
r49 = 'H2 + CO <=> HCO + H'
r50 = 'HO2 + CO <=> CO2 + OH'

Reaction{
   'HNO + H <=> H2 + NO',
   fr={'Arrhenius', A=0.48e13, n=0.0, C=0.0},
   label='r51'
}

Reaction{
   'HNO + OH <=> H2O + NO',
   fr={'Arrhenius', A=0.36e14, n=0.0, C=0.0},
   label='r52'
}

r53 = 'NO + CO <=> CO2 + N'
r54 = 'NO2 + CO <=> NO + CO2'

Reaction{
   'NO + HO2 <=> HNO + O2',
   fr={'Arrhenius', A=0.72e13, n=0.5, C=5500.0},
   label='r55'
}

Reaction{
   'HNO + O <=> NO + OH',
   fr={'Arrhenius', A=0.5e12, n=0.5, C=0.0},
   label='r56'
}

Reaction{
   'HNO3 + O <=> HO2 + NO2',
   fr={'Arrhenius', A=0.10e12, n=0.0, C=0.0},
   label='r57'
}

Reaction{
   'HO2 + NO2 <=> HNO2 + O2',
   fr={'Arrhenius', A=0.20e12, n=0.0, C=0.0},
   label='r58'
}

r59 = 'HCO + O2 <=> CO + HO2'

Reaction{
   'O3 + HO2 <=> 2 O2 + OH',
   fr={'Arrhenius', A=0.10e12, n=0.0, C=1409.0},
   label='r60'
}

if options.H2_O2_only then
   includedReactions = {
      'r1', 'r2', 'r3', 'r4', 'r9', 'r12',
      'r14', 'r15', 'r16', 'r17', 'r18', 'r19',
      'r20', 'r21', 'r22', 'r23', 'r24', 'r25',
      'r26', 'r27', 'r28', 'r29', 'r30', 'r31',
      'r32', 'r34', 'r35', 'r60' }
   selectOnlyReactionsWithLabel(includedReactions)
end


