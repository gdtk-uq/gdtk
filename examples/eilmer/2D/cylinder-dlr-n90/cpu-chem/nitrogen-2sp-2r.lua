-- nitrogen-2sp-2r.lua
--
-- This chemical kinetic system provides
-- a simple nitrogen dissociation mechanism.
--
-- Author: Rowan J. Gollan
-- Date: 13-Mar-2009 (Friday the 13th)
-- Place: NIA, Hampton, Virginia, USA
--
-- History:
--   24-Mar-2009 - reduced file to minimum input
--   11-Aug-2015 - updated for dlang module

Reaction{
   'N2 + N2 <=> N + N + N2',
   fr={'Arrhenius', A=7.0e21,  n=-1.6, C=113200.0},
   br={'Arrhenius', A=1.09e16, n=-0.5, C=0.0}
}

Reaction{
   'N2 + N <=> N + N + N',
   fr={'Arrhenius', A=3.0e22,  n=-1.6, C=113200.0},
   br={'Arrhenius', A=2.32e21, n=-1.5, C=0.0}
}

