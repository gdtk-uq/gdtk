/**
 * mass_diffusion.d
 *
 * This module houses the models for mass diffusion
 * that can be coupled to the flow solver.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2017-06-06, first cut
 *
 */
module mass_diffusion;

import std.math;
import std.stdio;

import gas.gas_model;
import gas.physical_constants;
import flowstate;
import flowgradients;
import fvcore;

private immutable double SMALL_DIFFUSION_COEFFICIENT = 1.0e-20;

enum MassDiffusionModel { none, ficks_first_law }

string massDiffusionModelName(MassDiffusionModel i)
{
    final switch (i) {
    case MassDiffusionModel.none: return "none";
    case MassDiffusionModel.ficks_first_law: return "ficks_first_law";
    }
}

MassDiffusionModel massDiffusionModelFromName(string name)
{
    switch (name) {
    case "none": return MassDiffusionModel.none;
    case "ficks_first_law": return MassDiffusionModel.ficks_first_law;
    default: return MassDiffusionModel.none;
    }
}

interface MassDiffusion {
    void update_mass_fluxes(const FlowState fs, const FlowGradients grad, double[] jx, double[] jy, double[] jz);
}

MassDiffusion initMassDiffusion(GasModel gmodel, MassDiffusionModel mass_diffusion_model)
{
    switch (mass_diffusion_model) {
    case MassDiffusionModel.ficks_first_law:
	return new FicksFirstLaw(gmodel, true);
    default:
	throw new FlowSolverException("Selected mass diffusion model is not available.");
    }
}

class FicksFirstLaw : MassDiffusion {
    this(GasModel gmodel, bool withMassFluxCorrection=true)
    {
	_withMassFluxCorrection = withMassFluxCorrection;
	_gmodel = gmodel;
	_nsp = gmodel.n_species;
	
	_sigma.length = _nsp;
	_eps.length = _nsp;
	_D.length = _nsp;
	_M.length = _nsp;
	foreach (isp; 0 .. _nsp) {
	    _sigma[isp].length = _nsp;
	    _eps[isp].length = _nsp;
	    _D[isp].length = _nsp;
	    _M[isp].length = _nsp;
	}
	_D_avg.length = _nsp;
	_molef.length = _nsp;

	// Compute M_ij terms
	foreach (isp; 0 .. _nsp) {
	    foreach (jsp; 0 .. _nsp) {
		if (isp == jsp) continue;
		_M[isp][jsp] = 1.0/gmodel.mol_masses[isp] + 1.0/gmodel.mol_masses[jsp];
		_M[isp][jsp] = 2.0/_M[isp][jsp];
		_M[isp][jsp] *= 1.0e3; // from kg/mol to g/mol
	    }
	}
	// Compute sigma_ij terms
	foreach (isp; 0 .. _nsp) {
	    foreach (jsp; 0 .. _nsp) {
		if (isp == jsp) continue;
		_sigma[isp][jsp] = 0.5*(gmodel.LJ_sigmas[isp] + gmodel.LJ_sigmas[jsp]);
	    }
	}
	// Compute eps_ij terms
	foreach (isp; 0 .. _nsp) {
	    foreach (jsp; 0 .. _nsp) {
		if (isp == jsp) continue;
		_eps[isp][jsp] = sqrt(gmodel.LJ_epsilons[isp] * gmodel.LJ_epsilons[jsp]);
	    }
	}
    }
    void update_mass_fluxes(const FlowState fs, const FlowGradients grad, double[] jx, double[] jy, double[] jz)
    {
	_gmodel.massf2molef(fs.gas, _molef);
	computeBinaryDiffCoeffs(fs.gas.Ttr, fs.gas.p);
	computeAvgDiffCoeffs();
	foreach (isp; 0 .. _nsp) {
	    jx[isp] = -fs.gas.rho * _D_avg[isp] * grad.massf[isp][0];
	    jy[isp] = -fs.gas.rho * _D_avg[isp] * grad.massf[isp][1];
	    jz[isp] = -fs.gas.rho * _D_avg[isp] * grad.massf[isp][2];
	}
	if (_withMassFluxCorrection) {
	    // Correction as suggested by Sutton and Gnoffo, 1998  
	    double sum_x = 0.0;
	    double sum_y = 0.0;
	    double sum_z = 0.0;
	    foreach (isp; 0 .. _nsp) {
		sum_x += jx[isp];
		sum_y += jy[isp];
		sum_z += jz[isp];
	    }
	    foreach (isp; 0 .. _nsp) {
		jx[isp] = jx[isp] - fs.gas.massf[isp] * sum_x;
		jy[isp] = jy[isp] - fs.gas.massf[isp] * sum_y;
		jz[isp] = jz[isp] - fs.gas.massf[isp] * sum_z;
	    }
	}
    }

private:
    GasModel _gmodel;
    size_t _nsp;
    bool _withMassFluxCorrection;
    double[][] _sigma;
    double[][] _eps;
    double[][] _D;
    double[][] _M;
    double[] _D_avg;
    double[] _molef;
    // coefficients for diffusion collision integral calculation
    double _a = 1.06036;
    double _b = 0.15610;
    double _c = 0.19300;
    double _d = 0.47635;
    double _e = 1.03587;
    double _f = 1.52996;
    double _g = 1.76474;
    double _h = 3.89411;

    void computeBinaryDiffCoeffs(double T, double p)
    {
	// Expression from:
	// Reid et al.
	// The Properties of Gases and Liquids
	foreach (isp; 0 .. _nsp) {
	    foreach (jsp; 0 .. _nsp) {
		if (isp == jsp) continue;
		double T_star = T/_eps[isp][jsp];
		double omega = _a/(pow(T_star, _b));
		omega += _c/(exp(_d*T_star));
		omega += _e/(exp(_f*T_star));
		omega += _g/(exp(_h*T_star));

		_D[isp][jsp] = 1.0/(p/P_atm)/sqrt(_M[isp][jsp])/(_sigma[isp][jsp]*_sigma[isp][jsp])/omega;
		_D[isp][jsp] *= 0.00266*sqrt(T*T*T);
		_D[isp][jsp] *= 1.0e-4; // cm^2/s --> m^2/s
	    }
	}
    }

    void computeAvgDiffCoeffs()
    {
	foreach (isp; 0 .. _nsp) {
	    double sum = 0.0;
	    foreach (jsp; 0 .. _nsp) {
		if (isp == jsp) continue;
		// The following two if-statements should generally catch the
		// same flow condition, namely, a zero or very small presence of
		// a certain species.  In this case the diffusion is effectively
		// zero and its contribution to the mixture diffusion coefficient
		// may be ignored.
		//
		// The two statements are used for extra security in detecting the
		// condition.
		if (_D[isp][jsp] < SMALL_DIFFUSION_COEFFICIENT ) continue;  // there is effectively nothing to diffuse
		if (_molef[jsp] < SMALL_MOLE_FRACTION ) continue; // there is effectively nothing to diffuse
		sum += _molef[jsp] / _D[isp][jsp];
	    }
	    if (sum <= 0.0) {
		_D_avg[isp] = 0.0;
	    }
	    else {
		_D_avg[isp] = (1.0 - _molef[isp])/sum;
	    }
	}
    }
}




