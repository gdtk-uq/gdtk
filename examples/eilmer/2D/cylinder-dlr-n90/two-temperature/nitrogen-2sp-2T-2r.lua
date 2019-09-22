-- nitrogen-2sp-2T-2r.lua
--
-- Author: Rowan J. Gollan
-- Date: 2019-09-22

Reaction{
   'N2 + N2 <=> N + N + N2',
   fr={'Park', A=7.0e21,  n=-1.6, C=113200.0, s=0.5},
   br={'Arrhenius', A=1.09e16, n=-0.5, C=0.0}
}

Reaction{
   'N2 + N <=> N + N + N',
   fr={'Park', A=3.0e22,  n=-1.6, C=113200.0, s=0.5},
   br={'Arrhenius', A=2.32e21, n=-1.5, C=0.0}
}

