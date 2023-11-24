/**
 * Authors: Rowan G. and Nick G.
 * Date: 2021-09-12
 *
 * rps = Reid, Prausnitz and Sherwood
 *
 * References:
 * Reid, Prausnitz and Sherwood (1977)
 * The Properties of Gases and Liquids, 3rd edition.
 * McGraw-Hill, New York.
 *
 * Neufeld, Janzen and Aziz (1972)
 * Empirical equations to calculate 16 of the transport collision integrals
 * $\Omega^{(l,s)*}$ for the Lennard-Jones (12-6) potential.
 * The Journal of Chemical Physics, 57(3), pp. 1100-1102.
 *
 */

module gas.diffusion.rps_diffusion_coefficients;

import std.math;
import std.conv : to;
import ntypes.complex;
import nm.number;


import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.binary_diffusion_coefficients;

class RPSDiffusionCoefficients : BinaryDiffusionCoefficients {
public:
    this(in double[] mol_masses, in double[] LJ_sigmas, in double[] LJ_epsilons)
    {
        mNspecies = mol_masses.length;
        mEpsilons.length = mNspecies;
        mSigmas.length = mNspecies;
        mM.length = mNspecies;
        foreach (isp; 0 .. mNspecies) {
            mEpsilons[isp].length = mNspecies;
            mSigmas[isp].length = mNspecies;
            mM[isp].length = mNspecies;
        }
        foreach (isp; 0 .. mNspecies) {
            foreach (jsp; isp+1 .. mNspecies) {
                mEpsilons[isp][jsp] = sqrt(LJ_epsilons[isp] * LJ_epsilons[jsp]);
                mSigmas[isp][jsp] = 0.5*(LJ_sigmas[isp] + LJ_sigmas[jsp]);
                mM[isp][jsp] = 2.0/(1.0/mol_masses[isp] + 1.0/mol_masses[jsp])*1.0e3;
            }
        }
    }

    void compute_bdc(ref const(GasState) Q, ref number[][] D)
    {
        // Parameters from Neufeld et al. expression.
        // See eqn (4.67) in RJG PhD thesis.
        immutable double A = 1.06036;
        immutable double B = 0.15610;
        immutable double C = 0.19300;
        immutable double d = 0.47635; // To avoid shadowing D
        immutable double E = 1.03587;
        immutable double F = 1.52996;
        immutable double G = 1.76474;
        immutable double H = 3.89411;

        number T = Q.T;
        number p = Q.p;
        foreach (isp; 0 .. mNspecies) {
            foreach (jsp; isp+1 .. mNspecies) {
                // Compute collision integral with Neufeld et al fit.
                // Eqn (4.67) in RJG Thesis
                number T_star = T/mEpsilons[isp][jsp];
                number Omega = A/(pow(T_star, B));
                Omega += C/(exp(d*T_star));
                Omega += E/(exp(F*T_star));
                Omega += G/(exp(H*T_star));

                number numer = 0.00266*sqrt(T*T*T);
                numer *= 1.0e-4; // cm^2/s --> m^2/s
                number denom = p/P_atm;
                denom *= sqrt(mM[isp][jsp]);
                denom *= mSigmas[isp][jsp] * mSigmas[isp][jsp];
                denom *= Omega;
                D[isp][jsp] = numer/denom;
                D[jsp][isp] = D[isp][jsp];
            }
        }
    }

private:
    size_t mNspecies;
    double[][] mEpsilons;
    double[][] mSigmas;
    double[][] mM;
}
