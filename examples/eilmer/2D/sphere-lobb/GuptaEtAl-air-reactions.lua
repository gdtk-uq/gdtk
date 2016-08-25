--
-- This chemical kinetic system provides
-- a full 20 reaction scheme for air.
--
-- Reference:
-- ---------
-- Gupta, R.N., Yos, J.M., Thompson, R.A. and Lee, K.-P. (1990)
-- A Review of Reaction Rates and Thermodynamic and Transport
-- Properties for an 11-Species Air Model for Chemical and
-- Thermal Nonequilibrium Calculations to 30 000 K.
-- NASA Reference Publication 1232
--
-- Author: Fabian Zander
-- Date: 12-Mar-2010
-- Place: Munich, Germany

-- Species are N, O, N2, O2, NO, N+, O+, N2+, O2+, NO+, e-
-- If NO_SPECIES == 5 then species are:
--   N2, O2, N, O, NO
-- If NO_SPECIES == 7 then species are:
--   N2, O2, N, O, NO, NO+, e-
-- 
-- History
-- -------
-- 2016-01-19 : Updated by RJG for use with Dlang project

NO_SPECIES = 5

-- Currently two schemes are available in this file using 
-- the WITH_IONIZATION flag. If 'false' then the first
-- 6 reactions are used. If 'true' the full 20 reactions
-- are considered. If other schemes are required this can
-- be configured in the same manner using the 
-- 'select_reactions_by_label' function as demonstrated at
-- the end of the file.

WITH_IONIZATION = true

if NO_SPECIES == 5 then
   WITH_IONIZATION = false
end

-- We know that some of the species in the efficiencies list do NOT
-- apply when we don't consider ionization.  We choose to suppress
-- those warnings.

SUPPRESS_WARNINGS = true

-- Temperature limits.
-- These temperature limits are used to limit the evaluation of
-- the rate coefficients. Gupta et al give a vailidity range
-- up to 30,000K, so it may not be wise to evaluate these
-- rate coefficients for higher temperatures. Also, the
-- rate values do get meaningless at very low temperatures.
-- For this reason, we suggest the limits we set here.
-- WARNING: We have experienced a flow solver/chemistry coupling
-- issue for strong expansions when the lower limit is set too low.
-- 
Config{
   tempLimits = {lower=300.0, upper=30000.0}
}

Reaction{
	'O2 + M <=> O + O + M',
	fr={'Arrhenius', A=3.610e+18, n=-1.00, C=59400.00},
	br={'Arrhenius', A=3.010e+15, n=-0.50, C=0.0},
	label='r1',
	efficiencies={O2=9.0, N2=2.0, O=25.0, N=1.0, NO=1.0, 
		['NO+']=0.0, ['O2+']=0.0, ['N2+']=0.0, ['O+']=0.0, ['N+']=0.0, ['e-']=0.0}
}

-- The efficiency of N is set to 0.0 so that it does NOT
-- participate in reaction r2.
-- The dissociation of N2 due to collisions with N has
-- a different reaction rate, and is included as the
-- subsequent (r3) reaction entry.

Reaction{
	'N2 + M <=> N + N + M',
	fr={'Arrhenius', A=1.920e+17, n=-0.50, C=113100.00},
	br={'Arrhenius', A=1.090e+16, n=-0.50, C=0.0},
	label='r2',
	efficiencies={O2=1.0, N2=2.5, O=1.0, N=0.0, NO=1.0,
		['NO+']=0.0, ['O2+']=0.0, ['N2+']=0.0, ['O+']=0.0, ['N+']=0.0, ['e-']=0.0}
}

Reaction{
	'N2 + N <=> N + N + N',
	fr={'Arrhenius', A=4.150e+22, n=-1.50, C=113100.00},
	br={'Arrhenius', A=2.320e+21, n=-1.50, C=0.0},
	label='r3'
}

Reaction{
	'NO + M <=> N + O + M',
	fr={'Arrhenius', A=3.970e+20, n=-1.50, C=75600.00},
	br={'Arrhenius', A=1.010e+20, n=-1.50, C=0.0},
	label='r4',
	efficiencies={O2=1.0, N2=1.0, O=20.0, N=20.0, NO=20.0,
		['NO+']=0.0, ['O2+']=0.0, ['N2+']=0.0, ['O+']=0.0, ['N+']=0.0, ['e-']=0.0}
}

Reaction{
	'NO + O <=> O2 + N',
	fr={'Arrhenius', A=3.180e+09, n=1.00, C=19700.00},
	br={'Arrhenius', A=9.630e+11, n=0.50, C=3600.0},
	label='r5'
}

Reaction{
	'N2 + O <=> NO + N',
	fr={'Arrhenius', A=6.750e+13, n=0.00, C=37500.00},
	br={'Arrhenius', A=1.500e+13, n=0.00, C=0.0},
	label='r6'
}

Reaction{
	'N + O <=> NO+ + e-',
	fr={'Arrhenius', A=9.030e+09, n=0.50, C=32400.00},
	br={'Arrhenius', A=1.800e+19, n=-1.00, C=0.0},
	label='r7'
}

Reaction{
	'O + e- <=> O+ + e- + e-',
	fr={'Arrhenius', A=3.600e+31, n=-2.91, C=158000.00},
	br={'Arrhenius', A=2.200e+40, n=-4.50, C=0.0},
	label='r8'
}

Reaction{
	'N + e- <=> N+ + e- + e-',
	fr={'Arrhenius', A=1.100e+32, n=-3.14, C=169000.00},
	br={'Arrhenius', A=2.200e+40, n=-4.50, C=0.0},
	label='r9'
}

Reaction{
	'O + O <=> O2+ + e-',
	fr={'Arrhenius', A=1.600e+17, n=-0.98, C=80800.00},
	br={'Arrhenius', A=8.020e+21, n=-1.50, C=0.0},
	label='r10'
}

Reaction{
	'O + O2+ <=> O2 + O+',
	fr={'Arrhenius', A=2.920e+18, n=-1.11, C=28000.00},
	br={'Arrhenius', A=7.800e+11, n=0.50, C=0.0},
	label='r11'
}

Reaction{
	'N2 + N+ <=> N + N2+',
	fr={'Arrhenius', A=2.020e+11, n=0.81, C=13000.00},
	br={'Arrhenius', A=7.800e+11, n=0.50, C=0.0},
	label='r12'
}

Reaction{
	'N + N <=> N2+ + e-',
	fr={'Arrhenius', A=1.400e+13, n=0.00, C=67800.00},
	br={'Arrhenius', A=1.500e+22, n=-1.50, C=0.0},
	label='r13'
}

Reaction{
	'O2 + N2 <=> NO + NO+ + e-',
	fr={'Arrhenius', A=1.380e+20, n=-1.84, C=141000.00},
	br={'Arrhenius', A=1.000e+24, n=-2.50, C=0.0},
	label='r14'
}

Reaction{
	'NO + M <=> NO+ + e- + M',
	fr={'Arrhenius', A=2.200e+15, n=-0.35, C=100800.00},
	br={'Arrhenius', A=2.200e+26, n=-2.50, C=0.0},
	label='r15',
	efficiencies={O2=4.0, N2=1.0, O=0.0, N=0.0, NO=0.0,
		['NO+']=0.0, ['O2+']=0.0, ['N2+']=0.0, ['O+']=0.0, ['N+']=0.0, ['e-']=0.0}
}

Reaction{
	'O + NO+ <=> NO + O+',
	fr={'Arrhenius', A=3.630e+15, n=-0.60, C=50800.00},
	br={'Arrhenius', A=1.500e+13, n=0.00, C=0.0},
	label='r16'
}

Reaction{
	'N2 + O+ <=> O + N2+',
	fr={'Arrhenius', A=3.400e+19, n=-2.0, C=23000.00},
	br={'Arrhenius', A=2.480e+19, n=-2.20, C=0.0},
	label='r17'
}

Reaction{
	'N + NO+ <=> NO + N+',
	fr={'Arrhenius', A=1.000e+19, n=-0.93, C=61000.00},
	br={'Arrhenius', A=4.800e+14, n=0.00, C=0.0},
	label='r18'
}

Reaction{
	'O2 + NO+ <=> NO + O2+',
	fr={'Arrhenius', A=1.800e+15, n=0.17, C=33000.00},
	br={'Arrhenius', A=1.800e+13, n=0.50, C=0.0},
	label='r19'
}

Reaction{
	'O + NO+ <=> O2 + N+',
	fr={'Arrhenius', A=1.340e+13, n=0.31, C=77270.00},
	br={'Arrhenius', A=1.000e+14, n=0.00, C=0.0},
	label='r20'
}

if NO_SPECIES == 5 then
   selectOnlyReactionsWithLabel({'r1', 'r2', 'r3', 'r4', 'r5', 'r6'})
end

if NO_SPECIES == 7 then
   selectOnlyReactionsWithLabel({'r1', 'r2', 'r3', 'r4', 'r5', 'r6',
			      'r7', 'r14', 'r15'})
end
