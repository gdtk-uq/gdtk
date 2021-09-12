/**
 * Authors: Rowan G. and Nick G.
 * Date: 2021-09-12
 */

module gas.diffusion.binary_diffusion_coefficients;

import gas.gas_state;
import nm.complex;
import nm.number;

interface BinaryDiffusionCoefficients {
    @nogc void compute_bdc(GasState Q, ref number[][] D);
}

