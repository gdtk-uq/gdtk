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

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;


interface ExchangeChemistryCoupling {
    ExchangeChemistryCoupling dup();
    @nogc number Gvanish(in GasState gs);
    @nogc number Gappear(in GasState gs);
}

class ImpartialDissociation : ExchangeChemistryCoupling {
    /*
        Simplified coupling model that assumes dissociation is equally likely from any vibrational
        level. Essentially it's the model proposed in Treanor and Marrone, 1962, but written with 
        the nomenclature of Knab et al. 1995
    */
    this(lua_State *L) {
        this.D = getDouble(L, -1, "D");
        this.Thetav = getDouble(L, -1, "Thetav");
    }

    this(double D, double Thetav) {
        this.D = D;
        this.Thetav = Thetav;
    }

    ImpartialDissociation dup() {
        return new ImpartialDissociation(D, Thetav);
    }
    
    @nogc
    number Gvanish(in GasState gs) {
    /*
        Equation (38) from Knab, 1995, with alpha=1 and A=D, and U=-inf which greatly simplifies the expression.
        Note their explanation of these assumptions on the paragraph just above the equation, in the right
        hand column of text.
    */  
        number T = gs.T;
        number Tv = gs.T_modes[0];
        number Lambda = 1.0/(1.0/Tv - 1.0/T);
        return L(Lambda, D);
    }

    @nogc
    number Gappear(in GasState gs) {
    /*
        This function is derived from Gvanish as Tv -> T. Computed using a Taylor series expansion
        about 1/Lambda = 0. See Treanor and Marrone equation 13 for the idea.
    */  
        double Gapp = 0.5*(D - Thetav*R_universal);
        return to!number(Gapp);
    }

private:
    const double D, Thetav;

    @nogc const number L(number T, double Y){
    /*
        Equation (39) from Knab, 1995, the molar vibrational energy content of the harmonic oscillator.
    */
        return R_universal*Thetav/(exp(Thetav/T) - 1.0) - Y/(exp(Y/R_universal/T) - 1.0);
    }
}

ExchangeChemistryCoupling createExchangeChemistryCoupling(lua_State *L)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "ImpartialDissociation":
        return new ImpartialDissociation(L);
    default:
        string msg = format("The exchange chemistry coupling model: %s is not known.", model);
        throw new Error(msg);
    }
}

