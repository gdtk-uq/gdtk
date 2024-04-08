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

Reaction{
   'H2 + O2 <=> OH + OH',
   fr={"Arrhenius", A=0.170e14, n=0, C=24230.0},
   label='r1'
}

Reaction{
   'O2 + H <=> OH + O',
   fr={"Arrhenius", A=0.142e15, n=0, C=8250.0},
   label='r2'
}

Reaction{
   'H2 + OH <=> H2O + H',
   fr={"Arrhenius", A=0.316e08, n=1.8, C=1525.0},
   label='r3'
}

Reaction{
   'H2 + O <=> OH + H',
   fr={"Arrhenius", A=0.207e15, n=0, C=6920.0},
   label='r4'
}

Reaction{
   'OH + OH <=> H2O + O',
   fr={"Arrhenius", A=0.550e14, n=0, C=3520.0},
   label='r5'
}

Reaction{
   'OH + H + M <=> H2O + M',
   fr={"Arrhenius", A=0.221e23, n=-2, C=0.0},
   label='r6'
}

Reaction{
   'H + H + M <=> H2 + M',
   fr={"Arrhenius", A=0.655e18, n=-1, C=0.0},
   label='r7'
}
