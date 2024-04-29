/*
Definitions for chemistry/vibration coupling models

References:
 - "Effect of Dissociation on the Rate of Vibrational Relaxation"
   Charles E. Treanor and Paul V. Marrone, The Physics of Fluids, Volume 5, Number 9, September 1962

 - "Chemical Relaxation with Preferential Dissociation from Excited Vibrational Levels"
   Paul V. Marrone and Charles E. Treanor, The Physics of Fluids, Volume 6, Number 9, September 1963

 - "Theory and Validation of a Physically Consistent Couplied Vibration-Chemistry-Vibration Model"
   O. Knab and H. H. Fruhauf and E. W. Masserschmid, Journal of Thermophysics and Heat Transfer,
   Volume 9, Number 2, April-June 1995

@author: Nick Gibbons (21/05/14)
*/

module kinetics.exchange_chemistry_coupling;

import std.string;
import std.math;
import std.conv;
import std.stdio;

import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;

immutable number iGammaSwitchThreshold = to!number(1e-7); // Determined by trial and error, see notes 26/05/21
immutable double iGammaOverflowThreshold = 0.001; // See NNG notes 01/09/22

interface ExchangeChemistryCoupling {
    ExchangeChemistryCoupling dup();
    @nogc number Gvanish(in GasState gs);
    @nogc number Gappear(in GasState gs);
}

@nogc number L(number iT, double Y, double Thetav) {
    /*
    Equation (39) from Knab, 1995, the molar vibrational energy content of the
    harmonic oscillator.

    Notes: We use the inverse pseudo-temperature iT=1/T to save some division operations
    */
    return R_universal*Thetav/(exp(Thetav*iT) - 1.0) - Y/(exp(Y/R_universal*iT) - 1.0);
}

@nogc number L2(number iT, double Y, double Thetav){
    /*
    Equation (39) from Knab, 1995, the molar vibrational energy content of the harmonic oscillator.

    Notes: This expression is a limit as Y/Ru*iT -> large, which catches a NaN that otherwise
           appears due to exp(Y/R_u*iT) overflowing. It's technically possible, but very
           unlikely, that exp(Thetav*iT) can overflow as well, but we have no check for that
           presently. (NNG 01/09/22)
    */
    return R_universal*Thetav/(exp(Thetav*iT) - 1.0) - Y*exp(-Y/R_universal*iT);
}

@nogc number LTaylorSeries1(number iT, double Y, double Thetav){
    /*
    First order Taylor series expansion of equation 39, to handle precision 
    loss issues near iT=0.
    */
    return 0.5*(Y - Thetav*R_universal) + iT*(Thetav*Thetav*R_universal*R_universal - Y*Y)/12.0/R_universal;
}

class ImpartialDissociation : ExchangeChemistryCoupling {
    /*
        Simplified coupling model that assumes dissociation is equally likely from any vibrational
        level. Essentially it's the model proposed in Treanor and Marrone, 1962, but written with 
        the nomenclature of Knab et al. 1995
    */
    this(lua_State *L, int mode) {
        this.D = getDouble(L, -1, "D");
        this.Thetav = getDouble(L, -1, "Thetav");
        this.mode = mode;
    }

    this(double D, double Thetav, int mode) {
        this.D = D;
        this.Thetav = Thetav;
        this.mode = mode;
    }

    ImpartialDissociation dup() {
        return new ImpartialDissociation(D, Thetav, mode);
    }
    
    @nogc
    number Gvanish(in GasState gs) {
    /*
        Equation (38) from Knab, 1995, with alpha=1 and A=D, and U=-inf which greatly simplifies the expression.
        Note their explanation of these assumptions on the paragraph just above the equation, in the right
        hand column of text.
    */  
        number T = gs.T;
        number Tv = gs.T_modes[mode];
        number iGamma = (T - Tv)/(T*Tv); // Inverse of pseudo-temperature in equation (36)

        // Equation (39) gives garbage results when T and Tv are close to one another, so we switch
        // to a non-singular approximation if iGamma gets too small.
        if (fabs(iGamma)<iGammaSwitchThreshold){
            return LTaylorSeries1(iGamma, D, Thetav);
        // Equation (39) can also overflow if the translational temp is much higher than the Tv.
        // L2 has been derived using a limit that turns the troublesome exp(x) into an exp(-x),
        // which safely underflows to zero instead of a nasty floating point infinity.
        } else if (iGamma>iGammaOverflowThreshold) {
            return L2(iGamma, D, Thetav);
        } else {
            return L(iGamma, D, Thetav);
        }
    }

    @nogc
    number Gappear(in GasState gs) {
    /*
        This function is derived from Gvanish as Tv -> T. Computed using a Taylor series expansion
        about 1/Gamma = 0. See Treanor and Marrone equation 13 for the idea.
    */  
        double L0 = 0.5*(D - Thetav*R_universal);
        return to!number(L0);
    }

private:
    const double D, Thetav;
    int mode;
}

class MarroneTreanorDissociation : ExchangeChemistryCoupling {
    this (double theta_v, double D, double U, int mode) {
        this._theta_v = theta_v;
        this._D = D;
        this._U = U;
        this._mode = mode;
    }

    this (lua_State *L, int mode) {
        double theta_v = getDouble(L, -1, "Thetav");
        double D = getDouble(L, -1, "D");
        double U = getDouble(L, -1, "U");
        this(theta_v, D, U, mode);
    }

    ExchangeChemistryCoupling dup() {
        return new MarroneTreanorDissociation(_theta_v, _D, _U, _mode);
    }

    @nogc
    number Gvanish(in GasState gs) {
        number T = gs.T;
        number Tv = gs.T_modes[_mode];
        number iGamma = 1./Tv - 1./T - 1./_U;

        // see notes for ImpartialDissociation
        if (fabs(iGamma) < iGammaSwitchThreshold) {
            return LTaylorSeries1(iGamma, _D, _theta_v);
        }
        else if (iGamma > iGammaOverflowThreshold) {
            return L2(iGamma, _D, _theta_v);
        }
        return L(iGamma, _D, _theta_v);
    }

    @nogc
    number Gappear(in GasState gs) {
        // compute the energy per dissociation, but with T_v = T_tr
        number T_F_inv = -1./_U;
        return R_universal * (_theta_v / (exp(_theta_v * T_F_inv) - 1.) - _D / (exp(_D / R_universal * T_F_inv) - 1.));
    }

private:
    int _mode;
    double _U, _theta_v, _D;
}

class ModifiedMarroneTreanorDissociation : ExchangeChemistryCoupling {
    this (double theta_v, double T_D, double aU, double Ustar, int mode, double non_boltzmann) {
        this._theta_v = theta_v;
        this._T_D = T_D;
        this._aU = aU;
        this._Ustar = Ustar;
        this._mode = mode;
        this._non_boltzmann = non_boltzmann;
    }

    this (lua_State *L, int mode) {
        double theta_v = getDouble(L, -1, "Thetav");
        double T_D = getDouble(L, -1, "T_D");
        double aU = getDouble(L, -1, "aU");
        double Ustar = getDouble(L, -1, "Ustar");
        double non_boltzmann = getDoubleWithDefault(L, -1, "non_boltzmann", 0.0);
        this(theta_v, T_D, aU, Ustar, mode, non_boltzmann);
    }

    ExchangeChemistryCoupling dup() {
        return new ModifiedMarroneTreanorDissociation(_theta_v, _T_D, _aU, _Ustar, _mode, _non_boltzmann);
    }

    @nogc
    number Gvanish(in GasState gs) {
        number T = gs.T;
        number Tv = gs.T_modes[_mode];
        number Uinv = _aU / T + 1 / _Ustar;
        number T_F_inv = 1./Tv - 1./T - Uinv;

        number boltzmann = R_universal * ( _theta_v / (exp(_theta_v * T_F_inv) - 1.) - _T_D / (exp(_T_D * T_F_inv) - 1.));
        return boltzmann + _non_boltzmann * R_universal * _T_D;
    }

    @nogc
    number Gappear(in GasState gs) {
        // compute the energy per dissociation, but with T_v = T_tr
        number T_F_inv = - _aU / gs.T - 1 / _Ustar;
        return R_universal * (_theta_v / (exp(_theta_v * T_F_inv) - 1.) - _T_D / (exp(_T_D * T_F_inv) - 1.));
    }

private:
    int _mode;
    double _aU, _Ustar, _theta_v, _T_D;
    double _non_boltzmann;
}


class ImpartialChem : ExchangeChemistryCoupling {
    /* 
     *  This exchange mechanism should allow a species to be born or destroyed
     *  without changing the temperature of the species in the absence of
     *  other exchange mechanisms.
     */
public:
    this (GasModel gmodel, int mode, int isp){
        this._isp = isp;
        this._gmodel = gmodel;
        this._mode = mode;
        this._gs_Tr = GasState(gmodel);
    }

    ExchangeChemistryCoupling dup() {
        return new ImpartialChem(_gmodel, _isp, _mode);
    }

    @nogc number Gvanish(in GasState gs){
        return _gmodel.energyPerSpeciesInMode(gs, _isp, _mode) * _gmodel.mol_masses[_isp];
    }

    @nogc number Gappear(in GasState gs){
        _gs_Tr.copy_values_from(gs);
        _gs_Tr.T_modes[_mode] = gs.T_modes[_mode];
        return _gmodel.energyPerSpeciesInMode(_gs_Tr, _isp, _mode) * _gmodel.mol_masses[_isp];
    }

private:
    GasModel _gmodel;
    GasState _gs_Tr;
    int _isp;
    int _mode;
}

ExchangeChemistryCoupling createExchangeChemistryCoupling(lua_State *L, GasModel gmodel, int mode, int isp)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "ImpartialDissociation":
        return new ImpartialDissociation(L, mode);
    case "MarroneTreanorDissociation":
        return new MarroneTreanorDissociation(L, mode);
    case "ModifiedMarroneTreanorDissociation":
        return new ModifiedMarroneTreanorDissociation(L, mode);
    case "ImpartialChem":
        return new ImpartialChem(gmodel, mode, isp);
    default:
        string msg = format("The exchange chemistry coupling model: %s is not known.", model);
        throw new Error(msg);
    }
}

