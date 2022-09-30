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
import gas.thermo.energy_modes;
import gas.physical_constants;


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
            return LTaylorSeries1(iGamma, D);
        // Equation (39) can also overflow if the translational temp is much higher than the Tv.
        // L2 has been derived using a limit that turns the troublesome exp(x) into an exp(-x),
        // which safely underflows to zero instead of a nasty floating point infinity.
        } else if (iGamma>iGammaOverflowThreshold) {
            return L2(iGamma, D);
        } else {
            return L(iGamma, D);
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
    immutable number iGammaSwitchThreshold = to!number(1e-7); // Determined by trial and error, see notes 26/05/21
    int mode;
    immutable double iGammaOverflowThreshold = 0.001; // See NNG notes 01/09/22

    @nogc const number L(number iT, double Y){
    /*
        Equation (39) from Knab, 1995, the molar vibrational energy content of the harmonic oscillator.

        Notes: We use the inverse pseudo-temperature iT=1/T to save some division operations
    */
        return R_universal*Thetav/(exp(Thetav*iT) - 1.0) - Y/(exp(Y/R_universal*iT) - 1.0);
    }

    @nogc const number L2(number iT, double Y){
    /*
        Equation (39) from Knab, 1995, the molar vibrational energy content of the harmonic oscillator.

        Notes: This expression is a limit as Y/Ru*iT -> large, which catches a NaN that otherwise
               appears due to exp(Y/R_u*iT) overflowing. It's technically possible, but very
               unlikely, that exp(Thetav*iT) can overflow as well, but we have no check for that
               presently. (NNG 01/09/22)
    */
        return R_universal*Thetav/(exp(Thetav*iT) - 1.0) - Y*exp(-Y/R_universal*iT);
    }

    @nogc const number LTaylorSeries1(number iT, double Y){
    /*
        First order Taylor series expansion of equation 39, to handle precision loss issues near iT=0.
    */
        return 0.5*(Y - Thetav*R_universal) + iT*(Thetav*Thetav*R_universal*R_universal - Y*Y)/12.0/R_universal;
    }
}

class ImpartialChemElectronElectronic : ExchangeChemistryCoupling {
public:
    this (GasModel gmodel, int isp, int mode){
        this.isp = isp;
        this.gmodel = gmodel;
        this.mode = mode;
    }

    ExchangeChemistryCoupling dup() {
        return new ImpartialChemElectronElectronic(gmodel, isp, mode);
    }

    @nogc number Gvanish(in GasState gs){
        // Species vanish with electron/electronic energy equal to the average
        // electron/electronic energy at the electron/electronic temperature.
        // Multiply by molar mass to convert J/kg -> J/mol
        return gmodel.energyPerSpeciesInMode(gs, isp, mode) * gmodel.mol_masses[isp];
    }

    @nogc number Gappear(in GasState gs){
        // Species appear with electron/electronic energy equal to the average
        // electron/electronic energy at the electron/electronic temperature.
        // Multiply by molar mass to convert J/kg -> J/mol
        return gmodel.energyPerSpeciesInMode(gs, isp, mode) * gmodel.mol_masses[isp];
    }

private:
    GasModel gmodel;
    int isp;
    int mode;
}

class ElectronImpactIonisation : ExchangeChemistryCoupling {
public:
    this (number ionisation_energy){
        // convert eV to J/mol
        this.ionisation_energy = ionisation_energy * electron_volt_energy * Avogadro_number;
    }

    this (lua_State* L){
        double ionisation_energy = getDouble(L, -1, "ionisation_energy");
        this(to!number(ionisation_energy));
    }

    ElectronImpactIonisation dup() {
        return new ElectronImpactIonisation(ionisation_energy);
    }

    @nogc number Gappear(in GasState gs){
        /// electrons appear with the energy of ionisation for the species that ionised
        return -ionisation_energy;
    }

    @nogc number Gvanish(in GasState gs){
        /// electrons vanish with the energy of ionisation for the species absorbing it
        return -ionisation_energy;
    }

private:
    number ionisation_energy;
}

ExchangeChemistryCoupling createExchangeChemistryCoupling(lua_State *L, int mode, GasModel gmodel, int isp)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "ImpartialDissociation":
        return new ImpartialDissociation(L, mode);
    case "ImpartialChemElectronElectronic":
        return new ImpartialChemElectronElectronic(gmodel, isp, mode);
    case "ElectronAtTranslationTemp":
        return new ImpartialChemElectronElectronic(gmodel, isp, -1);
    case "ElectronImpactIonisation":
        return new ElectronImpactIonisation(L);
    default:
        string msg = format("The exchange chemistry coupling model: %s is not known.", model);
        throw new Error(msg);
    }
}

