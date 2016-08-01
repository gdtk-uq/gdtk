-- Author: Rowan J. Gollan
-- Date: 02-Feb-2010
-- Place: Poquoson, Virginia, USA
--
-- Adapted from Python file: evans_schexnayder.py
--
-- This file provides four chemical kinetic descriptions
-- of hydrogen combustion.  You can select between the various
-- options below by setting the 'model' variable below to one of
-- the strings listed below.
--
-- REDUCED  : a 7-species, 8-reactions description of hydrogen
--            combustion in pure oxygen
-- PURE_O2  : a 7-species, 16-reactions description of hydrogen
--            combustion in pure oxygen
-- IN_AIR   : a 12-species, 25-reactions description of hydrogen
--            combustion in air (N2 and O2)
-- INERT_N2 : an 8-species, 16-reactions description of hydrogen
--            combustion in air with inert N2 (acting as diluent only).
--
-- The numbering of reactions in this file corresponds to
-- Table 1 in Evans and Schexnayder (1980).
--
-- Reference:
--  Evans, J.S. and Shexnayder Jr, C.J. (1980)
--  Influence of Chemical Kinetics and Unmixedness
--  on Burning in Supersonic Hydrogen Flames
--  AIAA Journal 18:2 pp 188--193
--
-- History:
--  07-Mar-2006 -- first prepared
--
-- Species used in REDUCED: O, O2, H, H2, H2O, OH, N2
-- Species used in PURE_O2: O, O2, H, H2, H2O, OH, HO2
-- Species used in IN_AIR: O, O2, N, N2, H, H2, H2O, HO2, OH, NO, NO2, HNO2
-- Species used in INERT_N2: O, O2, H, H2, H2O, OH, HO2, N2

options = {
   REDUCED=true,
   PURE_O2=true,
   IN_AIR=true,
   INERT_N2=true
}

-- User selects model here
model = 'INERT_N2'

-- Check that selection is valid
if options[model] == nil then
   print("User selected model: ", model)
   print("is not valid.")
   print("Valid models are:")
   for m,_ in pairs(options) do
      print(m)
   end
end

Config{
   odeStep = {method='alpha-qss'},
--   tightTempCoupling = true
}

Reaction{
   'HNO2 + M <=> NO + OH + M',
   fr={'Arrhenius', A=5.0e17, n=-1.0, C=25000.0},
   br={'Arrhenius', A=8.0e15, n=0.0, C=-1000.0},
   label='r1'
}

Reaction{
   'NO2 + M <=> NO + O + M',
   fr={'Arrhenius', A=1.1e16, n=0.0, C=32712.0},
   br={'Arrhenius', A=1.1e15, n=0.0, C=-941.0},
   label='r2'
}

Reaction{
   'H2 + M <=> H + H + M',
   fr={'Arrhenius', A=5.5e18, n=-1.0, C=51987.0},
   br={'Arrhenius', A=1.8e18, n=-1.0, C=0.0},
   label='r3'
}

Reaction{
   'O2 + M <=> O + O + M',
   fr={'Arrhenius', A=7.2e18, n=-1.0, C=59340.0},
   br={'Arrhenius', A=4.0e17, n=-1.0, C=0.0},
   label='r4'
}

Reaction{
   'H2O + M <=> OH + H + M',
   fr={'Arrhenius', A=5.2e21, n=-1.5, C=59386.0},
   br={'Arrhenius', A=4.4e20, n=-1.5, C=0.0},
   label='r5'
}

Reaction{
   'OH + M <=> O + H + M',
   fr={'Arrhenius', A=8.5e18, n=-1.0, C=50830.0},
   br={'Arrhenius', A=7.1e18, n=-1.0, C=0.0},
   label='r6'
}

Reaction{
   'HO2 + M <=> H + O2 + M',
   fr={'Arrhenius', A=1.7e16, n=0.0, C=23100.0},
   br={'Arrhenius', A=1.1e16, n=0.0, C=-440.0},
   label='r7'
}

Reaction{
   'H2O + O <=> OH + OH',
   fr={'Arrhenius', A=5.8e13, n=0.0, C=9059.0},
   br={'Arrhenius', A=5.3e12, n=0.0, C=503.0},
   label='r8'
}

Reaction{
   'H2O + H <=> OH + H2',
   fr={'Arrhenius', A=8.4e13, n=0.0, C=10116.0},
   br={'Arrhenius', A=2.0e13, n=0.0, C=2600.0},
   label='r9'
}

Reaction{
   'O2 + H <=> OH + O',
   fr={'Arrhenius', A=2.2e14, n=0.0, C=8455.0},
   br={'Arrhenius', A=1.5e13, n=0.0, C=0.0},
   label='r10'
}

Reaction{
   'H2 + O <=> OH + H',
   fr={'Arrhenius', A=7.5e13, n=0.0, C=5586.0},
   br={'Arrhenius', A=3.0e13, n=0.0, C=4429.0},
   label='r11'
}

Reaction{
   'H2 + O2 <=> OH + OH',
   fr={'Arrhenius', A=1.7e13, n=0.0, C=24232.0},
   br={'Arrhenius', A=5.7e11, n=0.0, C=14922.0},
   label='r12'
}

Reaction{
   'H2 + O2 <=> H + HO2',
   fr={'Arrhenius', A=1.9e13, n=0.0, C=24100.0},
   br={'Arrhenius', A=1.3e13, n=0.0, C=0.0},
   label='r13'
}

Reaction{
   'OH + OH <=> H + HO2',
   fr={'Arrhenius', A=1.7e11, n=0.5, C=21137.0},
   br={'Arrhenius', A=6.0e13, n=0.0, C=0.0},
   label='r14'
}

Reaction{
   'H2O + O <=> H + HO2',
   fr={'Arrhenius', A=5.8e11, n=0.5, C=28686.0},
   br={'Arrhenius', A=3.0e13, n=0.0, C=0.0},
   label='r15'
}

Reaction{
   'OH + O2 <=> O + HO2',
   fr={'Arrhenius', A=3.7e11, n=0.64, C=27840.0},
   br={'Arrhenius', A=1.0e13, n=0.0, C=0.0},
   label='r16'
}

Reaction{
   'H2O + O2 <=> OH + HO2',
   fr={'Arrhenius', A=2.0e11, n=0.5, C=36296.0},
   br={'Arrhenius', A=1.2e13, n=0.0, C=0.0},
   label='r17'
}

Reaction{
   'H2O + OH <=> H2 + HO2',
   fr={'Arrhenius', A=1.2e12, n=0.21, C=39815.0},
   br={'Arrhenius', A=1.7e13, n=0.0, C=12582.0},
   label='r18'
}

Reaction{
   'O + N2 <=> N + NO',
   fr={'Arrhenius', A=5.0e13, n=0.0, C=37940.0},
   br={'Arrhenius', A=1.1e13, n=0.0, C=0.0},
   label='r19'
}

Reaction{
   'H + NO <=> N + OH',
   fr={'Arrhenius', A=1.7e14, n=0.0, C=24500.0},
   br={'Arrhenius', A=4.5e13, n=0.0, C=0.0},
   label='r20'
}

Reaction{
   'O + NO <=> N + O2',
   fr={'Arrhenius', A=2.4e11, n=0.5, C=19200.0},
   br={'Arrhenius', A=1.0e12, n=0.5, C=3120.0},
   label='r21'
}

Reaction{
   'NO + OH <=> H + NO2',
   fr={'Arrhenius', A=2.0e11, n=0.5, C=15500.0},
   br={'Arrhenius', A=3.5e14, n=0.0, C=740.0},
   label='r22'
}


Reaction{
   'NO + O2 <=> O + NO2',
   fr={'Arrhenius', A=1.0e12, n=0.0, C=22800.0},
   br={'Arrhenius', A=1.0e13, n=0.0, C=302.0},
   label='r23'
}

Reaction{
   'NO2 + H2 <=> H + HNO2',
   fr={'Arrhenius', A=2.4e13, n=0.0, C=14500.0},
   br={'Arrhenius', A=5.0e11, n=0.5, C=1500.0},
   label='r24'
}

Reaction{
   'NO2 + OH <=> NO + HO2',
   fr={'Arrhenius', A=1.0e11, n=0.5, C=6000.0},
   br={'Arrhenius', A=3.0e12, n=0.5, C=1200.0},
   label='r25'
}

reactions_list = {}

if model == 'REDUCED' then
   reactions_list = {'r3', 'r4', 'r5', 'r6', 'r8', 'r9', 'r10', 'r11'}
end

if model == 'PURE_O2' or model == 'INERT_N2' then
   reactions_list = {'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10',
		     'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18'}
end


if model ~= 'IN_AIR' then
   -- For all other models we select only a subset.
   selectOnlyReactionsWithLabel(reactions_list)
end
