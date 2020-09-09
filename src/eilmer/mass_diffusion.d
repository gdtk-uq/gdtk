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
import std.conv;
import std.string;

import nm.complex;
import nm.number;
import util.lua;
import gas;

import globalconfig;
import flowstate;
import flowgradients;
import fvcore;

private immutable double SMALL_DIFFUSION_COEFFICIENT = 1.0e-20;

enum MassDiffusionModel { none, ficks_first_law }

@nogc
string massDiffusionModelName(MassDiffusionModel i)
{
    final switch (i) {
    case MassDiffusionModel.none: return "none";
    case MassDiffusionModel.ficks_first_law: return "ficks_first_law";
    }
}

@nogc
MassDiffusionModel massDiffusionModelFromName(string name)
{
    switch (name) {
    case "none": return MassDiffusionModel.none;
    case "ficks_first_law": return MassDiffusionModel.ficks_first_law;
    default: return MassDiffusionModel.none;
    }
}

interface MassDiffusion {
    @nogc
    void update_mass_fluxes(FlowState fs, const FlowGradients grad,
                            number[] jx, number[] jy, number[] jz);
}

MassDiffusion initMassDiffusion(GasModel gmodel,
                                bool sticky_electrons,
                                MassDiffusionModel mass_diffusion_model,
                                bool withConstantLewisNumber,
                                double Lewis,
                                bool withSpeciesSpecificLewisNumbers)
{
    switch (mass_diffusion_model) {
    case MassDiffusionModel.ficks_first_law:
        return new FicksFirstLaw(gmodel, sticky_electrons, true, withConstantLewisNumber,
                                 Lewis, withSpeciesSpecificLewisNumbers);
    default:
        throw new FlowSolverException("Selected mass diffusion model is not available.");
    }
}

class FicksFirstLaw : MassDiffusion {
    this(GasModel gmodel,
         bool sticky_electrons,
         bool withMassFluxCorrection=true,
         bool withConstantLewisNumber=false,
         double Lewis=1.0,
         bool withSpeciesSpecificLewisNumbers=false)
    {
        _withMassFluxCorrection = withMassFluxCorrection;
        _withConstantLewisNumber = withConstantLewisNumber;
        _withSpeciesSpecificLewisNumbers = withSpeciesSpecificLewisNumbers;
        _Le = Lewis;
        
        _gmodel = gmodel;
        _nsp = (sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
        
        _LeS.length = _nsp;
        _D.length = _nsp;
        foreach (isp; 0 .. _nsp) {
            _D[isp].length = _nsp;
        }
        _D_avg.length = _nsp;
        _molef.length = _nsp;

        if (!withConstantLewisNumber) { 
            if (withSpeciesSpecificLewisNumbers) {
                foreach (isp; 0 .. _nsp) {
                    _LeS[isp] = gmodel.Le[isp];
                }
            }
        }
    }

    @nogc
    void update_mass_fluxes(FlowState fs, const FlowGradients grad,
                            number[] jx, number[] jy, number[] jz)
    {
        version(multi_species_gas) {
            _gmodel.massf2molef(fs.gas, _molef);
            if (_withConstantLewisNumber) {
                number Cp = _gmodel.Cp(fs.gas);
                number alpha = fs.gas.k/(fs.gas.rho*Cp);
                foreach (isp; 0 .. _nsp) _D_avg[isp] = alpha/_Le;
            }
            else {
                if (!_withSpeciesSpecificLewisNumbers)
                    _gmodel.binary_diffusion_coefficients(fs.gas, _D);
                computeAvgDiffCoeffs(fs.gas);
            }
            foreach (isp; 0 .. _nsp) {
                jx[isp] = -fs.gas.rho * _D_avg[isp] * grad.massf[isp][0];
                jy[isp] = -fs.gas.rho * _D_avg[isp] * grad.massf[isp][1];
                jz[isp] = -fs.gas.rho * _D_avg[isp] * grad.massf[isp][2];
            }
            if (_withMassFluxCorrection) {
                // Correction as suggested by Sutton and Gnoffo, 1998  
                number sum_x = 0.0;
                number sum_y = 0.0;
                number sum_z = 0.0;
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
        } else {
            // this function is gitless for single-species gas
        }
    } // end update_mass_fluxes()

private:
    GasModel _gmodel;
    size_t _nsp;
    bool _withMassFluxCorrection;
    bool _withConstantLewisNumber;
    bool _withSpeciesSpecificLewisNumbers;
    double _Le = 1.0;
    number[][] _D;
    double[] _LeS;
    number[] _D_avg;
    number[] _molef;

    @nogc
    void computeAvgDiffCoeffs(GasState Q)
    {
        if(_withSpeciesSpecificLewisNumbers) {
            _gmodel.update_trans_coeffs(Q);
            number Prandtl = _gmodel.Prandtl(Q);
            foreach (isp; 0 .. _nsp) {
                _D_avg[isp] = Q.mu / (Q.rho * Prandtl * _LeS[isp]); 
            }
        }
        else {
            foreach (isp; 0 .. _nsp) {
                number sum = 0.0;
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
}

// This lua wrapper is somewhat fragile.
// It works presently (2018-09-11) because there
// is only one type of mass diffusion model available.
// This will need a re-work if it has wider utility.

// TODO: Put this in the lua gas model
//extern(C) int luafn_computeBinaryDiffCoeffs(lua_State *L)
//{
//    auto gmodel = GlobalConfig.gmodel_master;
//    auto n_species = (GlobalConfig.sticky_electrons) ? gmodel.n_heavy : gmodel.n_species;
//    // Expect temperature, pressure, and a table to push values into.
//    auto T = to!number(luaL_checknumber(L, 1));
//    auto p = to!number(luaL_checknumber(L, 2));
//
//    auto mdmodel = cast(FicksFirstLaw) GlobalConfig.massDiffusion;
//    mdmodel.computeBinaryDiffCoeffs(T, p);
//
//    foreach (isp; 0 .. n_species) {
//        lua_rawgeti(L, 3, isp);
//        foreach (jsp; 0 .. n_species) {
//            if (isp != jsp) {
//                lua_pushnumber(L, mdmodel._D[isp][jsp]);
//            }
//            else {
//                lua_pushnumber(L, 0.0);
//            }
//            lua_rawseti(L, -2, jsp);
//        }
//        lua_pop(L, 1);
//    }
//    return 0;
//}
//

