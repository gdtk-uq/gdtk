/* 
Modular Turbulence Modelling Interface

@author: Nick N. Gibbons (n.gibbons@uq.edu.au)
*/
import std.conv;
import std.math;
import std.stdio;
import std.string;
import std.json;
import flowstate;
import flowgradients;
import gas;
import json_helper;
import nm.number;
import nm.complex;

/*
Abstract base class defines functions all turbulence models must have:

*/
class TurbulenceModelObject{
    this() {} // For some reason having this(); instead of this() {} didn't work??? 
    this (const JSONValue config){}
    this (TurbulenceModelObject other){}

    // Property functions to be overridden.
    @nogc abstract string modelName() const;
    @nogc abstract size_t nturb() const;
    @nogc abstract bool isTurbulent() const;

    // Methods to be overridden.
    abstract TurbulenceModelObject dup();
    @nogc abstract void compute_source_terms(const FlowState fs,const FlowGradients grad, const number dwall, ref number[] source) const;
    @nogc abstract number turbulent_viscosity(const FlowState fs, const FlowGradients grad, const number dwall) const;
    @nogc abstract number turbulent_conductivity(const FlowState fs, GasModel gm) const;
    @nogc abstract number turbulent_signal_frequency(const FlowState fs) const;
    @nogc abstract number turbulent_kinetic_energy(const FlowState fs) const;
    @nogc abstract string primitive_variable_name(size_t i) const;
    @nogc abstract number primitive_variable_defaults(size_t i) const;
    @nogc abstract number viscous_transport_coefficient(size_t i) const;

    // Common methods
    override string toString() const
    {
    char[] repr;
    repr ~= "Model Name: " ~ modelName ~ "=(\n";
    repr ~= "n_equations=\"" ~ to!string(nturb) ~"\"\n";
    repr ~= "isTurbulent=\"" ~ to!string(isTurbulent) ~"\"\n";
    repr ~= ")";
    return to!string(repr);
    }

} // end class TurbulenceModelObject

/*
Object representing turbulence model in laminar flow.
 Notes: This should give sensible answers such as mu_t=0.0 if asked

*/
class noTurbulenceModel : TurbulenceModelObject {
    this (const JSONValue config){
        writeln("In no turbulence model class constructor.");
    }
    this (noTurbulenceModel other){
        writeln("In no turbulence model dup class constructor.");
    }

    // This seems weird but is apparently idiomatic?
    @nogc final override string modelName() const {return "none";}
    @nogc final override size_t nturb() const {return 0;}
    @nogc final override bool isTurbulent() const {return false;}

    final override noTurbulenceModel dup() {
        return new noTurbulenceModel(this);
    }

    @nogc final override
    void compute_source_terms(const FlowState fs,const FlowGradients grad, const number dwall, 
                              ref number[] source) const {
        return; 
    }

    @nogc final override number turbulent_viscosity(const FlowState fs, const FlowGradients grad, const number dwall) const {
        number mu_t = 0.0;
        return mu_t;
    }
    @nogc final override number turbulent_conductivity(const FlowState fs, GasModel gm) const {
        number k_t = 0.0;
        return k_t;
    }

    @nogc final override number turbulent_signal_frequency(const FlowState fs) const {
        number turb_signal = 1e-99; // We never want this to be the limiting factor (1/s)
        return turb_signal;
    }

    @nogc final override number turbulent_kinetic_energy(const FlowState fs) const {
        number tke=0.0;
        return tke;
    }

    @nogc final override string primitive_variable_name(size_t i) const {
        return "";
    }

    @nogc final override number primitive_variable_defaults(size_t i) const {
        number def = 0.0;
        return def;
    }
    @nogc final override number viscous_transport_coefficient(size_t i) const {
        number mu_eff = 0.0;
        return mu_eff;
    }
}


/*
Object representating the wilcox k-omega turbulence model.
 Notes:
  - Under construction
  - Original implementation by Wilson Chan et al.
  - See "Suitability of the k-w turbulence model for scramjet flowfield simulations"
    W. Y. K. Chan and P. A. Jacobs and D. J. Mee
    Int. J. Numer. Meth. Fluids (2011)

*/
class kwTurbulenceModel : TurbulenceModelObject {
    this (const JSONValue config){
        //writeln("In k_omega turbulence model class constructor.");

        // What to do about default values? inherit from globalconfig??? Throw an error if not found?
        number Pr_t = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        bool axisymmetric = getJSONbool(config, "axisymmetric", false);
        int dimensions = getJSONint(config, "dimensions", 2);
        this(Pr_t, axisymmetric, dimensions);
    }

    this (kwTurbulenceModel other){
        //writeln("In kw turbulence model dup class constructor.");
        // This constructir can access other's private variables
        // because they are the same class
        this(other.Pr_t, other.axisymmetric, other.dimensions);
        return;
    }

    this (number Pr_t, bool axisymmetric, int dimensions) {
        //writeln("In kw turbulence model specific class constructor");
        this.Pr_t = Pr_t;
        this.axisymmetric = axisymmetric;
        this.dimensions = dimensions;
        writeln("Pr_t =", Pr_t);
        writeln("axisymmetric =", axisymmetric);
        writeln("dimensions =", dimensions);
    }
    @nogc final override string modelName() const {return "k_omega";}
    @nogc final override size_t nturb() const {return 2;}
    @nogc final override bool isTurbulent() const {return true;}

    final override kwTurbulenceModel dup() {
        writeln("Called kw dup method");
        return new kwTurbulenceModel(this);
    }

    @nogc final override
    void compute_source_terms(const FlowState fs,const FlowGradients grad, const number dwall,
                              ref number[] source) const {
        /*
        Compute k-omega source terms.
        
        Production and Dissipation expressions for turbulence kinetic energy
        and turbulence frequency (or pseudo-vorticity). Based on Wilcox's 2006 model.
        
        Jan 2007: Initial implementation (Jan-Pieter Nap, PJ)
        Dec 2008: Implementation of the 3D terms (W Chan)
        Jan 2011: Minor modification to allow for implicit updating of tke and omega (W Chan)
                  All "fs->tke" and "fs->omega" instances are replaced with tke and omega.
        Jul 2014: Port to D by PJ
        Feb 2020: Moved here from fvcell "k_omega_time_derivatives" function (NNG)
        */
        // FIXME: You'll need a source_term_limits method 
        return; 
    }

    @nogc final override number turbulent_viscosity(const FlowState fs, const FlowGradients grad, const number dwall) const {
        /*
        Calculate mu_t, the turbulence viscosity.

        Notes:
          - Copied from fvcell.d function turbulence_viscosity_k_omega

        */
        /*
        number S_bar_squared;
        double C_lim = 0.875;
        double beta_star = 0.09;

        number dudx = grad.vel[0][0];
        number dudy = grad.vel[0][1];
        number dvdx = grad.vel[1][0];
        number dvdy = grad.vel[1][1];
        if ( dimensions == 2 ) {
            // 2D cartesian or 2D axisymmetric
            if ( axisymmetric ) {
                // 2D axisymmetric
                //number v_over_y = fs.vel.y / pos[0].y;
                number v_over_y = fs.vel.y / dwall;
                S_bar_squared = dudx*dudx + dvdy*dvdy + v_over_y*v_over_y
                    - 4.0/9.0 * (dudx + dvdy + v_over_y)
                    * (dudx + dvdy + v_over_y)
                    + 0.5 * (dudy + dvdx) * (dudy + dvdx) ;
            } else {
                // 2D cartesian
                S_bar_squared = dudx*dudx + dvdy*dvdy
                    - 4.0/9.0 * (dudx + dvdy) * (dudx + dvdy)
                    + 0.5 * (dudy + dvdx) * (dudy + dvdx);
            }
        } else {
            number dudz = grad.vel[0][2];
            number dvdz = grad.vel[1][2];
            number dwdx = grad.vel[2][0];
            number dwdy = grad.vel[2][1];
            number dwdz = grad.vel[2][2];
            S_bar_squared =  dudx*dudx + dvdy*dvdy + dwdz*dwdz
                - 4.0/9.0*(dudx + dvdy + dwdz)*(dudx + dvdy + dwdz)
                + 0.5 * (dudy + dvdx) * (dudy + dvdx)
                + 0.5 * (dudz + dwdx) * (dudz + dwdx)
                + 0.5 * (dvdz + dwdy) * (dvdz + dwdy);
        }
        S_bar_squared = fmax(0.0, S_bar_squared);
        number tke = fs.turb[0];
        number omega = fs.turb[1];
        number omega_t = fmax(omega, C_lim*sqrt(2.0*S_bar_squared/beta_star));
        number mu_t = fs.gas.rho * tke / omega_t;
        */
        number mu_t = 0.0;
        return mu_t;
    }
    @nogc final override number turbulent_conductivity(const FlowState fs, GasModel gm) const {
        // Warning: Make sure fs.mu_t is up to date before calling this method
        //number k_t = gm.Cp(fs.gas) * fs.mu_t / Pr_t;
        number k_t = 0.0;
        return k_t;
    }
    @nogc final override number turbulent_signal_frequency(const FlowState fs) const {
        number turb_signal = fs.turb[1];
        return turb_signal;
    }

    @nogc final override number turbulent_kinetic_energy(const FlowState fs) const {
        number tke=fs.turb[0];
        return tke;
    }

    @nogc final override string primitive_variable_name(size_t i) const {
        return _varnames[i];
    }

    @nogc final override number primitive_variable_defaults(size_t i) const {
        return _vardefaults[i];
    }

    @nogc final override number viscous_transport_coefficient(size_t i) const {
        // For line 657 of fvinterface bit where kw have slightly different transport coefficients 
        number mu_eff = 0.0;
        return mu_eff;
    }

private:
    immutable number Pr_t;
    immutable bool axisymmetric;
    immutable int dimensions;
    immutable string[2] _varnames = ["tke", "omega"];
    immutable number[2] _vardefaults = [0.0, 1.0];
}


TurbulenceModelObject init_turbulence_model(const string turbulence_model_name, const JSONValue config)
/*
Interface for generating polymorphic turbulence models, similar to ../gas/init_gas_model.d

@author: Nick Gibbons
*/
{
    TurbulenceModelObject turbulence_model;
    switch (turbulence_model_name) {
    case "none":
        turbulence_model = new noTurbulenceModel(config);
        break;
    case "k_omega":
        turbulence_model = new kwTurbulenceModel(config);
        break;
    default:
        string errMsg = format("The turbulence model '%s' is not available.", turbulence_model_name);
        throw new Error(errMsg);
    }
    return turbulence_model;
} // end init_turbulence_model()

