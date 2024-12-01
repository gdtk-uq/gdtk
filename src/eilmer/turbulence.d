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
import util.json_helper;
import nm.number;
import ntypes.complex;
import globalconfig;
import geom;
import fvcell;
import fvinterface;

/*
Abstract base class defines functions all turbulence models must have:
 - Consider adding an immutable int MAX_TURBULENT_EQUATIONS to control compile time sizing.

*/
class TurbulenceModel{
    this() {}
    this (const JSONValue config){}
    this (TurbulenceModel other){}

    // Property functions to be overridden.
    @nogc abstract string modelName() const;
    @nogc abstract size_t nturb() const;
    @nogc abstract bool isTurbulent() const;
    @nogc abstract bool needs_dwall() const;

    // Methods to be overridden.
    abstract TurbulenceModel dup();
    @nogc abstract void source_terms(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar, const number dwall, const number L_min, const number L_max, ref number[] source) const;
    @nogc abstract number turbulent_viscosity(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar, const number dwall) const;
    @nogc abstract number turbulent_conductivity(ref const(FlowState) fs, GasModel gm) const;
    @nogc abstract number turbulent_signal_frequency(ref const(FlowState) fs) const;
    @nogc abstract number turbulent_kinetic_energy(ref const(FlowState) fs) const;
    @nogc abstract number[3] turbulent_kinetic_energy_transport(ref const(FlowState) fs, ref const(FlowGradients) grad) const;
    @nogc abstract string primitive_variable_name(size_t i) const;
    abstract int primitive_variable_index(string tqName) const;
    @nogc abstract number turb_limits(size_t i) const;
    @nogc abstract number viscous_transport_coeff(ref const(FlowState) fs, size_t i) const;
    @nogc abstract bool is_valid(ref const(FlowStateLimits) fsl, const number[] turb) const;
    @nogc abstract number tke_rhoturb_derivatives(ref const(FlowState) fs, size_t i) const;
    @nogc abstract void set_flowstate_at_wall(const int gtl, const FVInterface IFace, const FVCell cell, ref FlowState fs) const;

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

} // end class TurbulenceModel

/*
Object representing turbulence model in laminar flow.
 Notes: This should give sensible answers such as mu_t=0.0 if asked

*/
class noTurbulenceModel : TurbulenceModel {
    this (){}
    this (const JSONValue config){}    // Used in GlobalConfig constructor
    this (noTurbulenceModel other){}   // Used in dup, in localconfig constructor

    // This seems weird but is apparently idiomatic?
    @nogc final override string modelName() const {return "none";}
    @nogc final override size_t nturb() const {return 0;}
    @nogc final override bool isTurbulent() const {return false;}
    @nogc final override bool needs_dwall() const {return false;}

    final override noTurbulenceModel dup() {
        return new noTurbulenceModel(this);
    }

    @nogc final override
    void source_terms(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar,
                      const number dwall, const number L_min, const number L_max,
                      ref number[] source) const {
        return;
    }

    @nogc final override number turbulent_viscosity(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar, const number dwall) const {
        number mu_t = 0.0;
        return mu_t;
    }
    @nogc final override number turbulent_conductivity(ref const(FlowState) fs, GasModel gm) const {
        number k_t = 0.0;
        return k_t;
    }

    @nogc final override number turbulent_signal_frequency(ref const(FlowState) fs) const {
        number turb_signal = 1e-99; // We never want this to be the limiting factor (1/s)
        return turb_signal;
    }

    @nogc final override number turbulent_kinetic_energy(ref const(FlowState) fs) const {
        number tke=0.0;
        return tke;
    }

    @nogc final override
    number[3] turbulent_kinetic_energy_transport(ref const(FlowState) fs, ref const(FlowGradients) grad) const {
        number[3] qtke;
        qtke[0] = 0.0;
        qtke[1] = 0.0;
        qtke[2] = 0.0;
        return qtke;
    }

    @nogc final override string primitive_variable_name(size_t i) const {
        return "";
    }

    final override int primitive_variable_index(string tqName) const {
        return -99;
    }

    @nogc final override number turb_limits(size_t i) const {
        number def = 0.0;
        return def;
    }
    @nogc final override number viscous_transport_coeff(ref const(FlowState) fs, size_t i) const {
        number mu_eff = 0.0;
        return mu_eff;
    }
    @nogc final override bool is_valid(ref const(FlowStateLimits) fsl, const number[] turb) const {
        return true;
    }
    @nogc final override number tke_rhoturb_derivatives(ref const(FlowState) fs, size_t i) const {
        number dtke_drhoturb = 0.0;
        return dtke_drhoturb;
    }

    @nogc final override
    void set_flowstate_at_wall(const int gtl, const FVInterface IFace, const FVCell cell, ref FlowState fs) const {
        /*
            Set the interface value of each turbulent primitive,
            given the nearest cell above it inside the flow domain
        */
        return;
    }
}

version(turbulence){
/*
Object representating the wilcox k-omega turbulence model.
 Notes:
  - Under construction
  - Original implementation by Wilson Chan et al.
  - See "Suitability of the k-w turbulence model for scramjet flowfield simulations"
    W. Y. K. Chan and P. A. Jacobs and D. J. Mee
    Int. J. Numer. Meth. Fluids (2011)
  - Added the vorticity-based source term option. KAD 2022-11-09
    This modification was originally suggested in ref. [1] for airfoils,
    and suggested for blunt-body flows in ref. [2]. It is referred to as the
    Wilcox k-omega Two-Equation Models with Vorticity Source Term in ref. [3].
    references:
    [1] Improved Two-Equation k-omega Turbulence Models for Aerodynamic Flows, F. R. Menter, NASA TM 103975 (1992)
    [2] Turbulence on Blunt Bodies at Reentry Speeds, S. Sefidbakht, Ohio State University Thesis (2011)
    [3] https://turbmodels.larc.nasa.gov/wilcox.html
*/
class kwTurbulenceModel : TurbulenceModel {
    this (){
        number Pr_t            = GlobalConfig.turbulence_prandtl_number;
        bool axisymmetric      = GlobalConfig.axisymmetric;
        int dimensions         = GlobalConfig.dimensions;
        number max_mu_t_factor = GlobalConfig.max_mu_t_factor;
        bool use_vorticity_based_source_term = false;
        this(Pr_t, axisymmetric, dimensions, max_mu_t_factor, use_vorticity_based_source_term);
    }

    this (const JSONValue config, bool use_vorticity_based_source_term = false){
        // What to do about default values? inherit from globalconfig??? Throw an error if not found?
        number Pr_t            = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        bool axisymmetric      = getJSONbool(config, "axisymmetric", false);
        int dimensions         = getJSONint(config, "dimensions", 2);
        number max_mu_t_factor = getJSONdouble(config, "max_mu_t_factor", 300.0);
        this(Pr_t, axisymmetric, dimensions, max_mu_t_factor, use_vorticity_based_source_term);
    }

    this (kwTurbulenceModel other){
        this(other.Pr_t, other.axisymmetric, other.dimensions, other.max_mu_t_factor, other.use_vorticity_based_source_term);
        return;
    }

    this (number Pr_t, bool axisymmetric, int dimensions, number max_mu_t_factor, bool use_vorticity_based_source_term) {
        this.Pr_t = Pr_t;
        this.axisymmetric = axisymmetric;
        this.dimensions = dimensions;
        this.max_mu_t_factor = max_mu_t_factor ;
        this.use_vorticity_based_source_term = use_vorticity_based_source_term;
	_varindices["tke"] = 0;
	_varindices["omega"] = 1;
    }
    @nogc final override string modelName() const {return "k_omega";}
    @nogc final override size_t nturb() const {return 2;}
    @nogc final override bool isTurbulent() const {return true;}
    @nogc final override bool needs_dwall() const {return false;}

    final override kwTurbulenceModel dup() {
        return new kwTurbulenceModel(this);
    }

    @nogc final override
    void source_terms(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar,
                      const number dwall, const number L_min, const number L_max,
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
        double alpha = 0.52;
        double beta_0 = 0.0708;
        number beta;
        double beta_star = 0.09;
        number P_K, D_K, P_W, D_W;
        number cross_diff;
        double sigma_d = 0.0;
        number WWS, X_w, f_beta;

        number dudx = grad.vel[0][0];
        number dudy = grad.vel[0][1];
        number dvdx = grad.vel[1][0];
        number dvdy = grad.vel[1][1];
        number dtkedx = grad.turb[0][0];
        number dtkedy = grad.turb[0][1];
        number domegadx = grad.turb[1][0];
        number domegady = grad.turb[1][1];
        number tke = fs.turb[0];
        number omega = fs.turb[1];
        if ( dimensions == 2 ) {
            // 2D cartesian or 2D axisymmetric
            if ( axisymmetric ) {
                // 2D axisymmetric
                //number v_over_y = fs.vel.y / pos[0].y;
                number v_over_y = fs.vel.y / ybar;
                if (use_vorticity_based_source_term) {
                    // TODO: we need to check if there should be a v_over_y term here. KAD 2022-11-09
                    number omega2 = (dvdx-dudy)*(dvdx-dudy);
                    P_K = fs.mu_t*omega2 - 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy);
                } else {
                    // JP.Nap correction from 03-May-2007 (-v_over_y in parentheses)
                    // P_K -= 0.6667 * mu_t * v_over_y * (dudx+dvdy-v_over_y);
                    // Wilson Chan correction to JP Nap's version (13 Dec 2008)
                    P_K = 2.0 * fs.mu_t * (dudx*dudx + dvdy*dvdy)
                        + fs.mu_t * (dudy + dvdx) * (dudy + dvdx)
                        - 2.0/3.0 * fs.mu_t * (dudx + dvdy + v_over_y)
                        * (dudx + dvdy + v_over_y)
                        + 2.0 * fs.mu_t * (v_over_y) * (v_over_y)
                        - 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy + v_over_y);
                }
                WWS = 0.25 * (dvdx - dudy) * (dvdx - dudy) * v_over_y ;
            } else {
                // 2D cartesian
                if (use_vorticity_based_source_term) {
                    number omega2 = (dvdx-dudy)*(dvdx-dudy);
                    P_K = fs.mu_t*omega2 - 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy);
                } else {
                    P_K = 1.3333 * fs.mu_t * (dudx*dudx - dudx*dvdy + dvdy*dvdy)
                        + fs.mu_t * (dudy + dvdx) * (dudy + dvdx)
                        - 0.66667 * fs.gas.rho * tke * (dudx + dvdy);
                }
                WWS = 0.0 ;
            }
            cross_diff = dtkedx * domegadx + dtkedy * domegady ;
        } else {
            number dudz = grad.vel[0][2];
            number dvdz = grad.vel[1][2];
            number dwdx = grad.vel[2][0];
            number dwdy = grad.vel[2][1];
            number dwdz = grad.vel[2][2];
            number dtkedz = grad.turb[0][2];
            number domegadz = grad.turb[1][2];
            // 3D cartesian
            if (use_vorticity_based_source_term) {
                number omega2 = (dwdy-dvdz)*(dwdy-dvdz) + (dudz-dwdx)*(dudz-dwdx) + (dvdx-dudy)*(dvdx-dudy);
                P_K = fs.mu_t*omega2 - 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy + dwdz);
            } else {
                P_K = 2.0 * fs.mu_t * (dudx*dudx + dvdy*dvdy + dwdz*dwdz)
                    - 2.0/3.0 * fs.mu_t * (dudx + dvdy + dwdz) * (dudx + dvdy + dwdz)
                    - 2.0/3.0 * fs.gas.rho * tke * (dudx + dvdy + dwdz)
                    + fs.mu_t * (dudy + dvdx) * (dudy + dvdx)
                    + fs.mu_t * (dudz + dwdx) * (dudz + dwdx)
                    + fs.mu_t * (dvdz + dwdy) * (dvdz + dwdy) ;
            }
            cross_diff = dtkedx * domegadx + dtkedy * domegady + dtkedz * domegadz ;
            WWS = 0.25 * (dudy - dvdx) * (dudy - dvdx) * dwdz
                + 0.25 * (dudz - dwdx) * (dudz - dwdx) * dvdy
                + 0.25 * (dvdz - dwdy) * (dvdz - dwdy) * dudx
                + 0.25 * (dudy - dvdx) * (dvdz - dwdy) * (dwdx + dudz)
                + 0.25 * (dudz - dwdx) * (dwdy - dvdz) * (dudy + dvdx)
                + 0.25 * (dvdx - dudy) * (dudz - dwdx) * (dwdy + dvdx) ;
        } // end if myConfig.dimensions

        D_K = beta_star * fs.gas.rho * tke * omega;

        // Apply a limit to the tke production as suggested by Jeff White, November 2007.
        const double P_OVER_D_LIMIT = 25.0;
        P_K = fmin(P_K, P_OVER_D_LIMIT*D_K);

        if ( cross_diff > 0 ) sigma_d = 0.125;
        P_W = alpha * omega / fmax(tke, small_tke) * P_K +
            sigma_d * fs.gas.rho / fmax(omega, small_omega) * cross_diff;

        X_w = fabs(WWS / pow(beta_star*omega, 3)) ;
        f_beta = (1.0 + 85.0 * X_w) / (1.0 + 100.0 * X_w) ;
        beta = beta_0 * f_beta;
        D_W = beta * fs.gas.rho * omega * omega;

        source[0] = P_K - D_K;
        source[1] = P_W - D_W;

        // Note: This appears to be unused. Consider bringing it back if needed.
        // It will complicate the interface a bit...
        //if (myConfig.limit_tke_production) {
        //    // Apply a final limit on the rate of tke production, related to the local thermodynamic energy
        //    double dt = SimState.dt_global;
        //    double deltaT = myConfig.tke_production_limit_in_kelvins;
        //    number Cv = myConfig.gmodel.Cv(fs.gas);
        //    number maxRateBasedOnLimit = fs.gas.rho*Cv*deltaT/dt;
        //    Q_rtke = fmin(maxRateBasedOnLimit, Q_rtke);
        //}
        return;
    }

    @nogc final override number turbulent_viscosity(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar, const number dwall) const {
        /*
        Calculate mu_t, the turbulence viscosity.

        Notes:
          - Copied from fvcell.d function turbulence_viscosity_k_omega

        */
        double C_lim = 0.875;
        double beta_star = 0.09;
        number S_bar_squared;

        number dudx = grad.vel[0][0];
        number dudy = grad.vel[0][1];
        number dvdx = grad.vel[1][0];
        number dvdy = grad.vel[1][1];
        if ( dimensions == 2 ) {
            // 2D cartesian or 2D axisymmetric
            if ( axisymmetric ) {
                // 2D axisymmetric
                //number v_over_y = fs.vel.y / pos[0].y;
                number v_over_y = fs.vel.y / ybar;
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
        return mu_t;
    }
    @nogc final override number turbulent_conductivity(ref const(FlowState) fs, GasModel gm) const {
        // Warning: Make sure fs.mu_t is up to date before calling this method
        number k_t = gm.Cp(fs.gas) * fs.mu_t / Pr_t;
        return k_t;
    }
    @nogc final override number turbulent_signal_frequency(ref const(FlowState) fs) const {
        number turb_signal = fs.turb[1];
        return turb_signal;
    }

    @nogc final override number turbulent_kinetic_energy(ref const(FlowState) fs) const {
        number tke = fs.turb[0];
        return tke;
    }

    @nogc final override
    number[3] turbulent_kinetic_energy_transport(ref const(FlowState) fs, ref const(FlowGradients) grad) const {
        // k-w uses the same expression for \overline{rho uj" 0.5 ui" ui"} as in the tke equation
        // taken from fvinterface viscous_flux_calc
        number mu_effective = viscous_transport_coeff(fs, 0);
        mu_effective = fmin(mu_effective, max_mu_t_factor * fs.gas.mu);

        number[3] qtke;
        qtke[0] = mu_effective * grad.turb[0][0];
        qtke[1] = mu_effective * grad.turb[0][1];
        qtke[2] = mu_effective * grad.turb[0][2];
        return qtke;
    }

    @nogc final override string primitive_variable_name(size_t i) const {
        return _varnames[i];
    }

    final override int primitive_variable_index(string tqName) const {
        return _varindices.get(tqName, -99);
    }

    @nogc final override number turb_limits(size_t i) const {
        return _varlimits[i];
    }

    @nogc final override number viscous_transport_coeff(ref const(FlowState) fs, size_t i) const {
        /*
        For line 657 of fvinterface, where k and w have slightly different transport coefficients
        */
        number sigma = _sigmas[i];
        number mu_effective = fs.gas.mu + sigma * fs.gas.rho * fs.turb[0] / fs.turb[1];
        return mu_effective;
    }

    @nogc final override bool is_valid(ref const(FlowStateLimits) fsl, const number[] turb) const {
        bool isvalid = is_tke_valid(fsl, turb[0]) && is_omega_valid(turb[1]);
        return isvalid;
    }

    @nogc final override number tke_rhoturb_derivatives(ref const(FlowState) fs, size_t i) const {
        number dtke_drhoturb = _rhocoeffs[i];
        return dtke_drhoturb;
    }

    @nogc final override
    void set_flowstate_at_wall(const int gtl, const FVInterface IFace, const FVCell cell, ref FlowState fs) const {
        /*
            Set the interface value of each turbulent primitive,
            given the nearest cell above it inside the flow domain
        */
        if (cell.in_turbulent_zone) {
            number d0 = distance_between(cell.pos[gtl], IFace.pos);
            fs.turb[1] = ideal_omega_at_wall(cell, d0);
            fs.turb[0] = 0.0;
        } else {
            fs.turb[1] = cell.fs.turb[1];
            fs.turb[0] = cell.fs.turb[0];
        }
        return;
    }


private:
    immutable number Pr_t;
    immutable bool axisymmetric;
    immutable int dimensions;
    immutable number max_mu_t_factor;
    immutable string[2] _varnames = ["tke", "omega"];
    int[string] _varindices;
    immutable number[2] _varlimits = [0.0, 1.0];
    immutable number[2] _sigmas = [0.6, 0.5];
    immutable number[2] _rhocoeffs = [1.0, 0.0];
    immutable number small_tke = 0.1;
    immutable number small_omega = 1.0;
    immutable bool use_vorticity_based_source_term;

    @nogc bool is_tke_valid(ref const(FlowStateLimits) flowstate_limits, const number tke) const {
        if (!isFinite(tke.re)) {
            debug { writeln("Turbulence KE invalid number ", tke); }
            return false;
        }
        if (tke < flowstate_limits.min_tke) {
            debug { writeln("Turbulence KE below minimum ", tke); }
            return false;
        }
        if (tke > flowstate_limits.max_tke) {
            debug { writeln("Turbulence KE above maximum ", tke); }
            return false;
        }
        return true;
    }

    @nogc bool is_omega_valid(const number omega) const {
        if (!isFinite(omega.re)) {
            debug { writeln("Turbulence frequency invalid number ", omega); }
            return false;
        }
        if (omega <= 0.0) {
            debug { writeln("Turbulence frequency nonpositive ", omega); }
            return false;
        }
        return true;
    }

    @nogc const
    number ideal_omega_at_wall(const FVCell cell, number d0)
    // As recommended by Wilson Chan, we use Menter's correction
    // for omega values at the wall. This appears as Eqn A12 in
    // Menter's paper.
    // Reference:
    // Menter (1994)
    // Two-Equation Eddy-Viscosity Turbulence Models for
    // Engineering Applications.
    // AIAA Journal, 32:8, pp. 1598--1605
    // Notes: Moved here from boundary_interface_effect.d (NNG)
    {
        // Note: d0 is half_cell_width_at_wall.
        number nu = cell.fs.gas.mu / cell.fs.gas.rho;
        double beta1 = 0.075;
        return 10 * (6 * nu) / (beta1 * d0 * d0);
    }

}

/*
Object representing the Spalart Allmaras turbulence model
 Notes:
  - The model here is taken from "Modifications and Clarifications
    for the Implementation of the Spalart-Alllmaras Turbulence Model"
    by Steven R. Allmaras, Forrester T. Johnson and Philippe R. Spalart
    2012, ICCFD7 paper
 - Model constants and nomenclature from "Simulation and Dynamics of
   Hypersonic Turbulent Combustion", 2019 PhD Thesis by Nick Gibbons

 @author: Nick Gibbons (n.gibbons@uq.edu.au)
*/
class saTurbulenceModel : TurbulenceModel {
    this (){
        number Pr_t = GlobalConfig.turbulence_prandtl_number;
        this(Pr_t);
    }

    this (const JSONValue config){
        number Pr_t = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        this(Pr_t);
    }

    this (saTurbulenceModel other){
        this(other.Pr_t);
    }

    this (number Pr_t) {
        this.Pr_t = Pr_t;
	_varindices["nuhat"] = 0;
    }

    @nogc override string modelName() const {return "spalart_allmaras";}
    @nogc final override size_t nturb() const {return 1;}
    @nogc final override bool isTurbulent() const {return true;}
    @nogc override bool needs_dwall() const {return true;}

    override saTurbulenceModel dup() {
        return new saTurbulenceModel(this);
    }

    @nogc override
    void source_terms(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar,
                      const number dwall, const number L_min, const number L_max,
                      ref number[] source) const {
        /*
        Spalart Allmaras Source Terms:
        Notes:
           - This term consists of Production - Destruction +
             compressible dissipation, the blue, red, and magenta terms
             from Gibbons, 2019, equation 2.32
        */

        number nuhat = fs.turb[0];
        number rho = fs.gas.rho;
        number nu = fs.gas.mu/rho;
        number chi = nuhat/nu;
        number chi_cubed = chi*chi*chi;
        number fv1 = chi_cubed/(chi_cubed + cv1_cubed);
        number fv2 = 1.0 - chi/(1.0 + chi*fv1);
        number ft2 = compute_ft2(chi);
        number nut = nuhat*fv1;

        number d = compute_d(nut,nu,grad.vel,dwall,L_min,L_max,fv1,fv2,ft2);
        number Shat_by_nuhat = compute_Shat_mulitplied_by_nuhat(grad, nuhat, nu, d, fv1, fv2);
        number production = rho*cb1*(1.0 - ft2)*Shat_by_nuhat;

        number r = compute_r(Shat_by_nuhat, nuhat, d);
        number g = r + cw2*(pow(r,6.0) - r);
        number fw = (1.0 + cw3_to_the_sixth)/(pow(g,6.0) +  cw3_to_the_sixth);
        fw = g*pow(fw, 1.0/6.0);
        number destruction = rho*(cw1*fw - cb1/kappa/kappa*ft2)*(nuhat*nuhat/d/d);

        //// No axisymmetric corrections terms in dS/dxi dS/dxi
        number nuhat_gradient_squared = 0.0;
        foreach(i; 0 .. 3) nuhat_gradient_squared+=grad.turb[0][i]*grad.turb[0][i];
        number dissipation = cb2/sigma*rho*nuhat_gradient_squared;

        number T = production - destruction + dissipation;
        //debug {
        //    if (!isFinite(T.re)) {
        //        writeln("Turbulence source term is NaN.");
        //        writeln("  r=", r, " fw=", fw, " rho=", rho, " cw1=", cw1, " cb1=", cb1);
        //        writeln("  kappa=", kappa, " ft2=", ft2, " nuhat=", nuhat, " d=", d, " Shat_by_nuhat=", Shat_by_nuhat);
        //        writeln("  production=", production, " destruction=", destruction, " dissipation=", dissipation);
        //        writeln("  T=", T);
        //        throw new Error("Turbulence source term is not finite.");
        //    }
        //}
        source[0] = T;
        return;
    } // end source_terms()

    @nogc override number turbulent_viscosity(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar, const number dwall) const {
        /*
        Compute the turbulent viscosity mu_t from the SA transport variable nuhat
        See equation (1) from Allmaras (2012)
        */
        number nuhat = fs.turb[0];
        number rho= fs.gas.rho;
        number mu = fs.gas.mu;
        number chi = rho*nuhat/mu;
        number chi_cubed = chi*chi*chi;
        number fv1 = chi_cubed/(chi_cubed + cv1_cubed);
        number mu_t = rho*nuhat*fv1;
        return mu_t;
    }
    @nogc override number turbulent_conductivity(ref const(FlowState) fs, GasModel gm) const {
        // Warning: Make sure fs.mu_t is up to date before calling this method
        number k_t = gm.Cp(fs.gas) * fs.mu_t / Pr_t;
        return k_t;
    }

    @nogc final override number turbulent_signal_frequency(ref const(FlowState) fs) const {
        /*
        The SA model doesn't really have an equivalent frequency...
        Something to think about in the future
        */
        number turb_signal = 1e-99; // We never want this to be the limiting factor (1/s)
        return turb_signal;
    }

    @nogc final override number turbulent_kinetic_energy(ref const(FlowState) fs) const {
        /*
        SA model assumes tke is small, there is a way to estimate it if need be
        */
        number tke=0.0;
        return tke;
    }

    @nogc final override
    number[3] turbulent_kinetic_energy_transport(ref const(FlowState) fs, ref const(FlowGradients) grad) const {
        /*
        Since tke==0 everywhere...
        */
        number[3] qtke;
        qtke[0] = 0.0;
        qtke[1] = 0.0;
        qtke[2] = 0.0;
        return qtke;
    }

    @nogc final override string primitive_variable_name(size_t i) const {
        return _varnames[i];
    }

    final override int primitive_variable_index(string tqName) const {
        return _varindices.get(tqName, -99);
    }

    @nogc final override number turb_limits(size_t i) const {
        return _varlimits[i];
    }

    @nogc override number viscous_transport_coeff(ref const(FlowState) fs, size_t i) const {
        /*
        Viscous diffusion of nuhat, the green term in Gibbons (2019) equation 2.32
        */
        number nuhat = fs.turb[0];
        number rho = fs.gas.rho;
        number nu = fs.gas.mu/rho;
        number mu_eff = rho*(nu + nuhat)/sigma;
        return mu_eff;
    }

    @nogc final override bool is_valid(ref const(FlowStateLimits) fsl, const number[] turb) const {
        return is_nuhat_valid(turb[0]);
    }

    @nogc final override number tke_rhoturb_derivatives(ref const(FlowState) fs, size_t i) const {
        /*
        Since tke is zero everywhere...
        */
        number dtke_drhoturb = 0.0;
        return dtke_drhoturb;
    }

    @nogc final override
    void set_flowstate_at_wall(const int gtl, const FVInterface IFace, const FVCell cell, ref FlowState fs) const {
        /*
        nuhat is set to zero at no slip walls in accord with Allmaras (2012), eqn 7
        */
        fs.turb[0] = 0.0;
        fs.mu_t = 0.0;
        fs.k_t = 0.0;
        return;
    }

protected:
    immutable number Pr_t;
    immutable string[1] _varnames = ["nuhat"];
    int[string] _varindices;
    immutable number[1] _varlimits = [0.0];
    immutable double sigma = 2.0/3.0;
    immutable double cb1 = 0.1355;
    immutable double cb2 = 0.622;
    immutable double kappa = 0.41;
    immutable double cw1 = cb1/kappa/kappa + (1.0+cb2)/sigma;
    immutable double cw2 = 0.3;
    immutable double cw3 = 2.0;
    immutable double cw3_to_the_sixth = cw3*cw3*cw3*cw3*cw3*cw3;
    immutable double cv1 = 7.1;
    immutable double cv1_cubed = cv1*cv1*cv1;
    immutable double cv2 = 0.7;
    immutable double cv3 = 0.9;
    immutable double ct1 = 1.0;
    immutable double ct2 = 2.0;
    immutable double ct3 = 1.2;
    immutable double ct4 = 0.5;
    immutable double rlimit = 10.0;

    @nogc number
    compute_d(const number nut, const number nu, const number[3][3] velgrad,
              const number dwall, const number L_min, const number L_max,
              const number fv1,  const number fv2, const number ft2) const {
        return dwall;
    }

    @nogc number
    compute_ft2(const number chi) const pure {
        return ct3*exp(-ct4*chi*chi);
    }

    @nogc number
    compute_Shat_mulitplied_by_nuhat(ref const(FlowGradients) grad, const number nuhat, const number nu,
                 const number d, const number fv1, const number fv2) const pure {
        // No axisymmetric corrections since W is antisymmetric
        number Omega = compute_Omega(grad);

        // Clipping Equation: NNG 2.37, Allmaras (11/12)
        number Sbar = nuhat/kappa/kappa/d/d*fv2;
        number Shat;
        number Sthing;
        if (Sbar >= -cv2*Omega) {
            Shat = Omega + Sbar;
        } else {
            Sthing = Omega*(cv2*cv2*Omega + cv3*Sbar)/((cv3-2.0*cv2)*Omega - Sbar);
            Shat = Omega + Sthing;
        }
        return Shat*nuhat;
    }

    @nogc number
    compute_r(const number Shat_by_nuhat, const number nuhat, const number d) const pure {
    /*
        Notes: (21/06/18)
        In a complex numbers calculation, it is possible to get a very small value of
        Shat_by_nuhat that intermixes the real and complex components when calculating r.

        This can result in a very large negative r, which is is both unphysical and likely
        to overflow when computing fw. The correct behaviour is to set r = 10+0j as Shat
        gets small. Finite negative values of Shat should never happen, but even
        if they do we still want r = 10+0j, to keep fw in line. This could still malfunction
        for very small real values of nuhat, but if anyone tries that they will have larger
        problems than just this function.

        @author: Nick Gibbons
    */
        number l = nuhat*nuhat/kappa/kappa/d/d;
        number r;

        if (Shat_by_nuhat*rlimit<=l) {
            r = rlimit;
        } else {
            r = l/Shat_by_nuhat;
        }
        return r;
    }

    @nogc number
    compute_Omega(ref const(FlowGradients) grad) const pure{
        number Omega = 0.0;
        number Wij;
        foreach(i; 0 .. 3) {
            foreach(j; 0 .. 3) {
                 Wij = 0.5*(grad.vel[i][j] - grad.vel[j][i]);
                 Omega += Wij*Wij;
            }
        }
        Omega = sqrt(2.0*Omega);
        return Omega;
    }

    @nogc final bool is_nuhat_valid(const number nuhat) const {
        if (!isFinite(nuhat.re)) {
            debug { writeln("Turbulence nuhat invalid number ", nuhat); }
            return false;
        }
        if (nuhat < 0.0) {
            debug { writeln("Turbulence nuhat below minimum 0.0"); }
            return false;
        }
        return true;
    }
}
/*
    Spalart-Allmaras `BCM' variant for transitional flows:

    "A Revised One-Equation Transitional Model for External Aerodynamics",
    Cakmakcioglu, S. C., Bas, O., Mura, R., and Kaynak, U.
    AIAA Paper 2020-2706, June 2020, (10.2514/6.2020-2706)

    @author: Yu Chen and Nick Gibbons
*/
class sabcmTurbulenceModel : saTurbulenceModel {
    this (){
        number Pr_t = GlobalConfig.turbulence_prandtl_number;
        double Tu_inf = GlobalConfig.freestream_turbulent_intensity;
        this(Pr_t, Tu_inf);
    }

    this (const JSONValue config){
        number Pr_t = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        double Tu_inf = getJSONdouble(config, "freestream_turbulent_intensity", 0.01);
        this(Pr_t, Tu_inf);
    }

    this (sabcmTurbulenceModel other){
        this(other.Pr_t, other.Tu_inf);
    }

    this (number Pr_t, double Tu_inf) {
        this.Tu_inf = Tu_inf;
        super(Pr_t);
    }

    @nogc override string modelName() const {return "spalart_allmaras_bcm";}

    override sabcmTurbulenceModel dup() {
        return new sabcmTurbulenceModel(this);
    }

    @nogc override
    void source_terms(ref const(FlowState) fs, ref const(FlowGradients) grad, const number ybar,
                      const number dwall, const number L_min, const number L_max,
                      ref number[] source) const {
        /*
        Spalart Allmaras Source Terms:
        Notes:
           - SA production term modified by Yu Chen
             See: https://turbmodels.larc.nasa.gov/sa-bc_1eqn.html
        */

        number nuhat = fs.turb[0];
        number rho = fs.gas.rho;
        number nu = fs.gas.mu/rho;
        number chi = nuhat/nu;
        number chi_cubed = chi*chi*chi;
        number fv1 = chi_cubed/(chi_cubed + cv1_cubed);
        number fv2 = 1.0 - chi/(1.0 + chi*fv1);
        number ft2 = 0.0;  // no ft2 in sa-bcm
        number nut = nuhat*fv1;

        //additional parmeters for sa-bcm
        number mu = fs.gas.mu;
        number re_theta_c = 803.73*pow((Tu_inf*100.0+0.6067),-1.027);
        number chi1 = 0.002;
        number chi2 = 0.02;
        number Omega = compute_Omega(grad);

        number d = compute_d(nut,nu,grad.vel,dwall,L_min,L_max,fv1,fv2,ft2);
        number Shat_by_nuhat = compute_Shat_mulitplied_by_nuhat(grad, nuhat, nu, d, fv1, fv2);

        //additional parmeters for sa-bcm
        number re_nu = rho*d*d/mu*Omega;  //omega needs to be defined
        number re_theta = re_nu/2.193;
        number mu_t = turbulent_viscosity(fs, grad, ybar, dwall);
        number term1 = fmax(re_theta - re_theta_c, 0.0)/(chi1 * re_theta_c);
        number term2 = fmax(mu_t/(chi2*mu), 0.0);  //mu_t needs to be defined
        number gamma_bc = 1.0 - exp(-sqrt(term1)-sqrt(term2));
        number production = gamma_bc*rho*cb1*Shat_by_nuhat;  // Different terms to sa mdoel

        number r = compute_r(Shat_by_nuhat, nuhat, d);
        number g = r + cw2*(pow(r,6.0) - r);
        number fw = (1.0 + cw3_to_the_sixth)/(pow(g,6.0) +  cw3_to_the_sixth);
        fw = g*pow(fw, 1.0/6.0);
        number destruction = rho*cw1*fw*nuhat*nuhat/d/d;

        //// No axisymmetric corrections terms in dS/dxi dS/dxi
        number nuhat_gradient_squared = 0.0;
        foreach(i; 0 .. 3) nuhat_gradient_squared+=grad.turb[0][i]*grad.turb[0][i];
        number dissipation = cb2/sigma*rho*nuhat_gradient_squared;

        number T = production - destruction + dissipation;
        source[0] = T;
        return;
    } // end source_terms()
private:
    double Tu_inf; // freestream turbulence intensity
}

/*
Object representing the Spalart Allmaras turbulence model with constant wall distance

 @author: Daryl Bond (d.bond1@uq.edu.au)
*/

class sadTurbulenceModel : saTurbulenceModel {
    this (){
        number Pr_t = GlobalConfig.turbulence_prandtl_number;
        this(Pr_t);
    }

    this (const JSONValue config){
        number Pr_t = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        this(Pr_t);
    }

    this (sadTurbulenceModel other){
        this(other.Pr_t);
    }

    this (number Pr_t) {
        this.Pr_t = Pr_t;
    }

    @nogc override string modelName() const {return "spalart_allmaras_dwall_const";}

    override sadTurbulenceModel dup() {
        return new sadTurbulenceModel(this);
    }

protected:
    immutable number Pr_t;

    @nogc override number
    compute_d(const number nut, const number nu, const number[3][3] velgrad,
              const number dwall, const number L_min, const number L_max,
              const number fv1,  const number fv2, const number ft2) const {
        return to!number(1000);
    }
}


/*
Object representing the Modified Edwards version of the Spalart Allmaras model
 Notes:
 - The model here is taken from "Comparison of Eddy Viscosity-Transport Turbulence
   Models for Three-Dimensional, Shock-Separated Flowfields", Edwards and Chandra (1997)
 - This version is intended by a more stable in compressible flow, particularly for
   implicit solvers with large timesteps.
 - It also removes the ft2 sink term, which is designed to gently pull the nuhat
   values down in regions of no strain, to allow for tripping. The paper does not
   explain the reasoning for this, but I found leaving it in caused numerical problems.

 @author: Nick Gibbons (n.gibbons@uq.edu.au)
*/
class saeTurbulenceModel : saTurbulenceModel {
    this (){
        number Pr_t = GlobalConfig.turbulence_prandtl_number;
        this(Pr_t);
    }

    this (const JSONValue config){
        number Pr_t = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        this(Pr_t);
    }

    this (saeTurbulenceModel other){
        this(other.Pr_t);
    }

    this (number Pr_t) {
        super(Pr_t);
    }

    @nogc final override string modelName() const {return "spalart_allmaras_edwards";}

    final override saeTurbulenceModel dup() {
        return new saeTurbulenceModel(this);
    }

protected:
    @nogc override number
    compute_ft2(const number chi) const pure {
        number ft2 = 0.0;
        return ft2;
    }

    @nogc override number
    compute_Shat_mulitplied_by_nuhat(ref const(FlowGradients) grad, const number nuhat, const number nu,
                 const number d, const number fv1, const number fv2) const pure {
        /*
        The Edwards variant of Spalart-Allmaras contains a definition of Shat
        that is singular when nuhat is equal to zero, such as during the Newton Krylov
        startup phase where the boundary conditions are diffused into the flow.

        Since Shat*nuhat is always non-singular and both the Edwards and normal models
        need it, we compute that instead.
        */
        number S = 0.0;
        foreach(i; 0 .. 3) {
            foreach(j; 0 .. 3) {
                 S += (grad.vel[i][j] + grad.vel[j][i])*grad.vel[i][j];
            }
        }
        number div = grad.vel[0][0] + grad.vel[1][1] + grad.vel[2][2];
        S -= 2.0/3.0*div*div;
        number Shat_by_nuhat = sqrt(S)*(nu + fv1*nuhat);
        return Shat_by_nuhat;
    }

    @nogc override number
    compute_r(const number Shat_by_nuhat, const number nuhat, const number d) const pure {
        number normal_r = super.compute_r(Shat_by_nuhat, nuhat, d);
        number r = tanh(normal_r)/tanh(1.0);
        return r;
    }
}


/*
Object representing the IDDES high fidelity turbulence model
 Notes:
 - The model here is taken from "A Hybrid RANS-LES Approach with
   Delayed-DES and Wall Modelled LES Capabilities"
 - Model constants and nomenclature from "Simulation and Dynamics of
   Hypersonic Turbulent Combustion", 2019 PhD Thesis by Nick Gibbons

 @author: Nick Gibbons (n.gibbons@uq.edu.au)
*/
class iddesTurbulenceModel : saTurbulenceModel {
    this (){
        number Pr_t = GlobalConfig.turbulence_prandtl_number;
        this(Pr_t);
    }

    this (const JSONValue config){
        number Pr_t = getJSONdouble(config, "turbulence_prandtl_number", 0.89);
        this(Pr_t);
    }

    this (iddesTurbulenceModel other){
        this(other.Pr_t);
    }

    this (number Pr_t) {
        super(Pr_t);
    }

    @nogc final override string modelName() const {return "iddes";}

    final override iddesTurbulenceModel dup() {
        return new iddesTurbulenceModel(this);
    }
private:
    immutable double Cw = 0.15;
    immutable double CDES = 0.65;
    immutable double fwstar = 0.424;

protected:
    @nogc override number
    compute_d(const number nut, const number nu, const number[3][3] velgrad,
              const number dwall, const number L_min, const number L_max,
              const number fv1,  const number fv2, const number ft2) const {
    /*
        IDDES uses a much more complicated definition of the turbulent length scale d,
        that behaves like RANS close to walls and like LES far away, with some input
        from the flow itself to trigger the transition.

        See Equations 2.50 to 2.63 in Gibbons (2019)
    */
        number PSI_numerator = 1.0 - (cb1/cw1/kappa/kappa/fwstar*(ft2 + (1.0-ft2)*fv2));
        number PSI_denominator = fv1*fmax(1e-10, 1.0-ft2);
        number PSI = sqrt(fmin(1e2, PSI_numerator/PSI_denominator));

        number S = 0.0;
        foreach(i; 0 .. 3) foreach(j; 0 .. 3) S += velgrad[i][j]*velgrad[i][j];
        S = sqrt(S);
        number r_denominator = kappa*kappa*dwall*dwall*fmax(S, 1e-10);
        number rdt = nut/r_denominator;
        number rdl = nu/r_denominator;

        number ct2rdt = 1.63*1.63*rdt;
        number cl2rdl = 3.55*3.55*rdl;

        number ft = tanh(pow(ct2rdt, 3.0));
        number fl = tanh(pow(cl2rdl, 10.0));
        number fe2 = 1.0 - fmax(ft, fl);

        number alpha = 0.25 - dwall/L_max;
        number fe1 = (alpha >= 0.0) ? 2.0*exp(-11.09*alpha*alpha) : 2.0*exp(-9.0*alpha*alpha);
        number fe = fmax(fe1 - 1.0, 0.0)*PSI*fe2;

        number fB = fmin(2.0*exp(-9.0*alpha*alpha), 1.0);
        number fd = 1.0 - tanh(512.0*rdt*rdt*rdt);
        number fh = fmax(1.0 - fd, fB);

        number maxDelta = fmax(Cw*dwall, Cw*L_max);
        maxDelta = fmax(maxDelta, L_min);
        number dLES = fmin(maxDelta, L_max);
        number d = fh*(1.0 + fe)*dwall + (1.0 - fh)*CDES*PSI*dLES;
        return d;
    }
}


} // end version(turbulence)


TurbulenceModel init_turbulence_model(const string turbulence_model_name, const JSONValue config)
    /*
    Interface for generating polymorphic turbulence models, similar to ../gas/init_gas_model.d
       - JSONValue version

    @author: Nick Gibbons
    */
{
    TurbulenceModel turbulence_model;
    switch (turbulence_model_name) {
    case "none":
        turbulence_model = new noTurbulenceModel(config);
        break;
version(turbulence){
    case "k_omega":
        turbulence_model = new kwTurbulenceModel(config);
        break;
    case "k_omega_vorticity":
        turbulence_model = new kwTurbulenceModel(config, true);
        break;
    case "spalart_allmaras":
        turbulence_model = new saTurbulenceModel(config);
        break;
    case "spalart_allmaras_bcm":
        turbulence_model = new sabcmTurbulenceModel(config);
        break;
    case "spalart_allmaras_dwall_const":
        turbulence_model = new sadTurbulenceModel(config);
        break;
    case "spalart_allmaras_edwards":
        turbulence_model = new saeTurbulenceModel(config);
        break;
    case "spalart_allmaras_iddes":
        turbulence_model = new iddesTurbulenceModel(config);
        break;
}
    default:
        string errMsg = format("The turbulence model '%s' is not available.", turbulence_model_name);
        throw new Error(errMsg);
    }
    return turbulence_model;
} // end init_turbulence_model()

TurbulenceModel init_turbulence_model(const string turbulence_model_name)
    /*
    Interface for generating polymorphic turbulence models, similar to ../gas/init_gas_model.d
       - Default version. Be very careful with this, it just uses GlobalConfig values, whether
       or not they have been set correctly yet!

    @author: Nick Gibbons
    */
{
    TurbulenceModel turbulence_model;
    switch (turbulence_model_name) {
    case "none":
        turbulence_model = new noTurbulenceModel();
        break;
version(turbulence){
    case "k_omega":
        turbulence_model = new kwTurbulenceModel();
        break;
    case "k_omega_vorticity":
        turbulence_model = new kwTurbulenceModel();
        break;
    case "spalart_allmaras":
        turbulence_model = new saTurbulenceModel();
        break;
    case "spalart_allmaras_bcm":
        turbulence_model = new sabcmTurbulenceModel();
        break;
    case "spalart_allmaras_dwall_const":
        turbulence_model = new sadTurbulenceModel();
        break;
    case "spalart_allmaras_edwards":
        turbulence_model = new saeTurbulenceModel();
        break;
    case "spalart_allmaras_iddes":
        turbulence_model = new iddesTurbulenceModel();
        break;
}
    default:
        string errMsg = format("The turbulence model '%s' is not available.", turbulence_model_name);
        throw new Error(errMsg);
    }
    return turbulence_model;
} // end init_turbulence_model()

