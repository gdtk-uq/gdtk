/**
 * Authors: Rowan G. and Nick G.
 * Date: 2021-09-12
 */

module gas.diffusion.binary_diffusion_coefficients;

import gas.gas_state;
import ntypes.complex;
import nm.number;

interface BinaryDiffusionCoefficients {
    @nogc void compute_bdc(ref const(GasState) Q, ref number[][] D);
}

