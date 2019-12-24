-- Author: Wilson Chan, Fabs, Rowan, Peter J.
-- Date: 07-Oct-2010, 20-12-2019
-- Place: Centre for Hypersonics, UQ
--
-- References
-- ----------
-- Rogers, R.C. and Chinitz, W. (1983)
-- Using a Global Hydrogen-Air Combustion Model in Turbulent
-- Reacting Flow Calculations
-- AIAA Journal, 21:4, pp. 586--592
--
-- J Philip Drummond
-- A Two-Dimensional Numerical Simulation of a Supersonic
-- Chemically Reacting Mixing Layer
-- NASA Technical Memorandum 4055, December 1988.
-- 
-- Notes
-- -----
-- 1. The pre-exponential factor 'A' in the Rogers and Chinitz
-- rate is dependent on the equivalence ratio. Hence, the
-- global equivalence ratio should be set when using this model.
-- A value of phi=0.3 is used to match the Drummond diffuser calculation.
--
-- 2. The notation matches that in Phil Drummond's NASA TM.
-- The backward reaction rates have also come from that report.
-- Comment those br tables out to get the default backward rates.
--
-- 3. I (Rowan Gollan) don't recommend this scheme. It's troublesome
-- in a numerical sense. Rogers and Chinitz discuss some of this
-- in their paper.
--
-- 4. Mott's alpha-QSS is a good selection for stiff exothermic systems
-- but the default convergence tolerance needs tightening for this scheme.
Config{
   -- odeStep = {method='alpha-qss', eps1=5.0e-6, delta=1.0e-10, maxIters=10},
   odeStep = {method='rkf', errTol=1.0e-6},
   tightTempCoupling = true
}

phi = 0.3
A1 = (8.917*phi + 31.433/phi - 28.950)*1e47
A2 = (2.000 + 1.333/phi - 0.833*phi)*1e64

R_cal = 1.9872041 -- cal/(K.mol)

Reaction{
   'H2 + O2 <=> 2OH',
   fr={'Arrhenius', A=A1, n=-10, C=4865/R_cal},
   br={'Arrhenius', A=A1/26.164, n=-10, C=(4865/R_cal)-8992},
   label='r1'
}

Reaction{
   '2OH + H2 <=> 2H2O',
   fr={'Arrhenius', A=A2, n=-13, C=42500/R_cal},
   br={'Arrhenius', A=A2/2.682e-6, n=-14, C=(42500/R_cal)+69415},
   label='r2'
}


