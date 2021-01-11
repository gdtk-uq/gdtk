-- An input to prep-chem for dissociating nitrogen (two-temperature version)
-- Rates are from Kewley & Hornung (1974)
--
-- To prepare file for eilmer:
-- > prep-chem nitrogen-gas-model.lua KewleyHornung-nitrogen-dissociation.lua N2-2T-KH.chem
--
-- Reference:
-- Kewley & Hornung (1974)
-- Free-piston Shock-tube Study of Nitrogen Dissociation,
-- Chemical Physics Letters, 25(4), pp. 531--536
--
-- This file prepared by:
-- Rowan J. Gollan, 2021-01-07
--
-- NOTE: 'Park' refers to the means to compute the forward rate constant
--        in a 2-T context. The Park model is used which takes the rate-controlling
--        temperature in the Arrhenius expression as:
--        T_rate = (T_tr)^(s) * (T_v)^(1-s)
--
--        where s=0.5 here, giving a geometric average of temperatures,
--        but 's' may be varied.
--

Reaction{
   'N2 + N2 <=> N + N + N2',
   fr={'Park', A=2.3e29,  n=-3.5, C=113200.0, s=0.5}
}

Reaction{
   'N2 + N <=> N + N + N',
   fr={'Park', A=8.5e25,  n=-2.5, C=113200.0, s=0.5}
}
