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

import ntypes.complex;
import nm.number;
import util.lua;
import gas;

import globalconfig;
import flowstate;
import flowgradients;

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
    void update_mass_fluxes(ref FlowState fs, ref const(FlowGradients) grad,
                            number[] jx, number[] jy, number[] jz);
}

MassDiffusion initMassDiffusion(GasModel gmodel,
                                string diffusion_coefficient_type,
                                MassDiffusionModel mass_diffusion_model,
                                double Lewis)
{
    switch (mass_diffusion_model) {
    case MassDiffusionModel.ficks_first_law:
        return new FicksFirstLaw(gmodel, diffusion_coefficient_type, true,
                                 Lewis);
    default:
        throw new FlowSolverException("Selected mass diffusion model is not available.");
    }
}

class FicksFirstLaw : MassDiffusion {
    this(GasModel gmodel,
         string diffusion_coefficient_type,
         bool withMassFluxCorrection=true,
         double Lewis=1.0)
    {
        _withMassFluxCorrection = withMassFluxCorrection;
        diffusion_coefficient = initDiffusionCoefficient(gmodel, diffusion_coefficient_type, Lewis);
        //
        _gmodel = gmodel;
        _nsp = gmodel.n_species;
        _D_avg.length = gmodel.n_species;
    }

    @nogc
    void update_mass_fluxes(ref FlowState fs, ref const(FlowGradients) grad,
                            number[] jx, number[] jy, number[] jz)
    {
        version(multi_species_gas) {
            diffusion_coefficient.computeAvgDiffCoeffs(fs.gas, _gmodel, _D_avg);
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
        } // this function is gitless for single-species gas
    } // end update_mass_fluxes()

private:
    DiffusionCoefficient diffusion_coefficient;
    GasModel _gmodel;
    size_t _nsp;
    bool _withMassFluxCorrection;
    number[] _D_avg;
}

interface DiffusionCoefficient {
    @nogc
    void computeAvgDiffCoeffs(ref GasState Q, GasModel gm, ref number[] D_avg);
}

class ConstantLewisNumber : DiffusionCoefficient {
    this(size_t nsp, double Le) {
        this.Le = Le;
        this.nsp = nsp;
    }
    @nogc
    void computeAvgDiffCoeffs(ref const(GasState) Q, GasModel gmodel, ref number[] D_avg) {
        number Cp = gmodel.Cp(Q);
        number alpha = Q.k/(Q.rho*Cp);
        foreach (isp; 0 .. nsp) D_avg[isp] = alpha/Le;
    }
private:
    size_t nsp;
    double Le;
}

class SpeciesSpecificLewisNumbers : DiffusionCoefficient {
    this(size_t nsp, double[] Le) {
        if (Le.length==0) throw new Error("No Lewis numbers present in gas model");
        foreach(Le_isp; Le) this.LeS ~= Le_isp;
        this.nsp = nsp;
    }
    @nogc
    void computeAvgDiffCoeffs(ref const(GasState) Q, GasModel gmodel, ref number[] D_avg) {
        // [FIX-ME] gmodel.update_trans_coeffs(Q); // This feels bad. Shouldn't this be set already???
        number Prandtl = gmodel.Prandtl(Q);
        foreach (isp; 0 .. gmodel.n_species) {
            D_avg[isp] = Q.mu / (Q.rho * Prandtl * LeS[isp]);
        }
    }
private:
    size_t nsp;
    double[] LeS;
}

class BinaryDiffusion : DiffusionCoefficient {
    this(size_t nsp, bool is_plasma, double[] charge) {
        this.nsp = nsp;
        molef.length = nsp;
        D.length = nsp;
        foreach (isp; 0 .. nsp) D[isp].length = nsp;

        if (is_plasma) {
            with_ambipolar_diffusion = true;
            foreach (isp, icharge; charge){
                if (icharge>0.0) ion_idxs ~= isp;
                if (icharge<0.0) electron_idx = isp;
            }
        }
    }
    @nogc
    void computeAvgDiffCoeffs(ref const(GasState) Q, GasModel gmodel, ref number[] D_avg) {
        gmodel.massf2molef(Q, molef);
        gmodel.binary_diffusion_coefficients(Q, D);
        foreach (isp; 0 .. nsp) {
            number sum = 0.0;
            foreach (jsp; 0 .. nsp) {
                if (isp == jsp) continue;
                // The following two if-statements should generally catch the
                // same flow condition, namely, a zero or very small presence of
                // a certain species.  In this case the diffusion is effectively
                // zero and its contribution to the mixture diffusion coefficient
                // may be ignored.
                //
                // The two statements are used for extra security in detecting the
                // condition.
                if (D[isp][jsp] < SMALL_DIFFUSION_COEFFICIENT ) continue;  // there is effectively nothing to diffuse
                if (molef[jsp] < SMALL_MOLE_FRACTION ) continue; // there is effectively nothing to diffuse
                sum += molef[jsp] / D[isp][jsp];
            }
            if (sum <= 0.0) {
                D_avg[isp] = 0.0;
            }
            else {
                D_avg[isp] = (1.0 - molef[isp])/sum;
            }
        }
        if (with_ambipolar_diffusion) computeAmbipolarDiffusion(Q, D_avg);
        return;
    }

private:
    size_t nsp;
    number[] molef;
    number[][] D;
    size_t electron_idx;
    size_t[] ion_idxs;
    bool with_ambipolar_diffusion = false;
    immutable double epsilon = 1e-6;

    @nogc
    void computeAmbipolarDiffusion(ref const(GasState) Q, ref number[] D_avg) {
        /*
        Ambipolar diffusion is a correction to the diffusion process for charged particles that
        contrains electrons and ions to diffuse at the same rate, thus preserving charge neutrality.

        For details see chapter 5 of "Introduction to Plasma Physics and Controlled Fusion", Chen 2016.
        The expression here is a generalisation of equation his 5.18, rederived for a multispecies plasma
        by NNG. See notes from 21/03/11.
        */

        number DiZi = 0.0;
        number theta = 1.0;
        version(multi_T_gas){
            if (Q.T_modes.length>0) theta = Q.T_modes[$-1]/Q.T;
        }

        // The following modified definition of Z was designed to handle
        // numerical problems that occur when you have very tiny amounts of
        // ionisation.  See NNG notes from 22/02/01 for the derivation.
        foreach (isp; ion_idxs) {
            number Zi = (molef[isp]+epsilon/ion_idxs.length)/(molef[electron_idx]+epsilon);
            DiZi += D_avg[isp]*Zi;
        }
        number Da = DiZi*(1.0+theta)/(DiZi*theta/D_avg[electron_idx] + 1.0);

        foreach (isp; ion_idxs) {
            D_avg[isp] = Da;
        }
        D_avg[electron_idx] = Da;
        return;
    }
}

DiffusionCoefficient initDiffusionCoefficient(GasModel gmodel, string diffusion_coefficient_type, double Lewis)
{
    switch (diffusion_coefficient_type) {
    case "constant_lewis_number":
        return new ConstantLewisNumber(gmodel.n_species, Lewis);
    case "species_specific_lewis_numbers":
        return new SpeciesSpecificLewisNumbers(gmodel.n_species, gmodel.Le);
    case "binary_diffusion":
        return new BinaryDiffusion(gmodel.n_species, gmodel.is_plasma, gmodel.charge);
    case "none":
        throw new FlowSolverException("Diffusion model requires a valid diffusion_coefficient_type.");
    default:
        string errMsg=format("The diffusion_coefficient_type '%s' is not available.",diffusion_coefficient_type);
        throw new FlowSolverException(errMsg);
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
//    auto n_species = gmodel.n_species;
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

