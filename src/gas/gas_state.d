/**
 * The GasState class defines the data storage for a blob of gas.
 *
 * Authors: Rowan G. and Peter J.
 */

module gas.gas_state;

import std.stdio;
import std.conv;
import std.math : isFinite;
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.cea_gas;

struct GasState {
public:
    /// Thermodynamic properties.
    number rho;  /// density, kg/m**3
    number p;    /// presure, Pa
    number p_e;  /// electron pressure
    number a;    /// sound speed, m/s
    // For a gas in thermal equilibrium, all of the internal energies
    // are bundled together into u and are represented by a single
    // temperature T.
    number T;    /// temperature, K
    number u;    /// specific thermal energy, J/kg
    // For a gas in thermal nonequilibrium, the internal energies are
    // stored unbundled, with u being the trans-rotational thermal energy.
    // The array length will be determined by the specific model and,
    // to get the total internal energy,
    // the gas-dynamics part of the code will need to sum the array elements.
    number[] u_modes;  /// specific internal energies for thermal nonequilibrium model, J/kg
    number[] T_modes;  /// temperatures for internal energies for thermal nonequilibrium, K
    /// Transport properties
    number mu;   /// viscosity, Pa.s
    number k;  /// thermal conductivity for a single temperature gas, W/(m.K)
    number[] k_modes;  /// thermal conductivities for the nonequilibrium model, W/(m.K)
    number sigma;    /// electrical conductivity, S/m
    /// Composition
    number[] massf;  /// species mass fractions
    number[] rho_s;  /// species densities
    number quality;  /// vapour quality
    // A place to hang on to some CEA data, so that it can be called up
    // in CEAGas methods that don't have access to the original CEA output file.
    CEASavedData* ceaSavedData;

    @disable this();

    this(uint n_species, uint n_modes, bool includeSavedData=false)
    {
        massf.length     = n_species;
        rho_s.length     = n_species;
        u_modes.length   = n_modes;
        T_modes.length   = n_modes;
        k_modes.length   = n_modes;
        if (includeSavedData) { ceaSavedData = new CEASavedData; }
    }

    this(GasModel gm)
    {
        this(gm.n_species, gm.n_modes);
        if (cast(CEAGas)gm) { ceaSavedData = new CEASavedData; }
    }

    this(GasModel gm, in double p_init, in double T_init, in double[] T_modes_init,
         in double[] massf_init=[1.0,], in double quality_init=1.0,
         in double sigma_init=0.0)
    {
        p = p_init;
        p_e = p_init;
        T = T_init;
        size_t nsp = gm.n_species;
        size_t nmodes = gm.n_modes;
        T_modes.length = nmodes;
        foreach(i; 0 .. nmodes) { T_modes[i] = T_modes_init[i]; }
        u_modes.length = nmodes;
        k_modes.length = nmodes;
        massf.length = nsp;
        rho_s.length = nsp;
        foreach(i; 0 .. nsp) {
            massf[i] = (i < massf_init.length) ? massf_init[i] : 0.0;
        }
        quality = quality_init;
        sigma = sigma_init;
        if (cast(CEAGas)gm) { ceaSavedData = new CEASavedData; }
        // Now, evaluate the rest of the properties using the gas model.
        gm.update_thermo_from_pT(this);
        gm.update_sound_speed(this);
        gm.update_trans_coeffs(this);
        foreach(i; 0 .. nsp) {
            rho_s[i] = rho*massf[i];
        }
    }

    this(GasModel gm, in double p_init, in double T_init,
         in double[] massf_init=[1.0,], in double quality_init=1.0,
         in double sigma_init=0.0)
    {
        double[] T_modes;
        T_modes.length = gm.n_modes;
        foreach(ref Tmode; T_modes) { Tmode = T_init; }
        this(gm, p_init, T_init, T_modes, massf_init, quality_init, sigma_init);
    }

    this(in GasState other)
    {
        rho = other.rho;
        p = other.p;
        p_e = other.p_e;
        T = other.T;
        u = other.u;
        a = other.a;
        u_modes = other.u_modes.dup;
        T_modes = other.T_modes.dup;
        mu = other.mu;
        k = other.k;
        k_modes = other.k_modes.dup;
        sigma = other.sigma;
        massf = other.massf.dup;
        rho_s = other.rho_s.dup;
        quality = other.quality;
        if (other.ceaSavedData !is null) {
            ceaSavedData = new CEASavedData(*(other.ceaSavedData));
        }
    }

    @nogc void opAssign(in GasState other)
    {
        debug {
            if (u_modes.length != other.u_modes.length) { throw new Error("Incorrect u_modes length."); }
            if (T_modes.length != other.T_modes.length) { throw new Error("Incorrect T_modes length."); }
            if (k_modes.length != other.k_modes.length) { throw new Error("Incorrect k_modes length."); }
            if (massf.length != other.massf.length) { throw new Error("Incorrect massf length."); }
            if (rho_s.length != other.rho_s.length) { throw new Error("Incorrect rho_s length."); }
        }
        rho = other.rho;
        p = other.p;
        T = other.T;
        u = other.u;
        p_e = other.p_e;
        a = other.a;
        foreach (i; 0 .. u_modes.length) { u_modes[i] = other.u_modes[i]; }
        foreach (i; 0 .. T_modes.length) { T_modes[i] = other.T_modes[i]; }
        mu = other.mu;
        k = other.k;
        foreach (i; 0 .. k_modes.length) { k_modes[i] = other.k_modes[i]; }
        sigma = other.sigma;
        foreach (i; 0 .. massf.length) { massf[i] = other.massf[i]; }
        foreach (i; 0 .. rho_s.length) { rho_s[i] = other.rho_s[i]; }
        quality = other.quality;
    }

    @nogc void copy_values_from(ref const(GasState) other)
    {
        this = other;
    }

    @nogc void copy_average_values_from(ref const(GasState) gs0, ref const(GasState) gs1, double w0=0.5)
    // Avoids memory allocation, it's all in place.
    {
        double w1 = 1.0 - w0;
        rho = w0*gs0.rho + w1*gs1.rho;
        p = w0*gs0.p + w1*gs1.p;
        T = w0*gs0.T + w1*gs1.T;
        u = w0*gs0.u + w1*gs1.u;
        p_e = w0*gs0.p_e + w1*gs1.p_e;
        a = w0*gs0.a + w1*gs1.a;
        foreach(i; 0 .. u_modes.length) { u_modes[i] = w0*gs0.u_modes[i] + w1*gs1.u_modes[i]; }
        foreach(i; 0 .. T_modes.length) { T_modes[i] = w0*gs0.T_modes[i] + w1*gs1.T_modes[i]; }
        mu = w0*gs0.mu + w1*gs1.mu;
        k = w0*gs0.k + w1*gs1.k;
        foreach(i; 0 .. k_modes.length) { k_modes[i] = w0*gs0.k_modes[i] + w1*gs1.k_modes[i]; }
        sigma = w0*gs0.sigma + w1*gs1.sigma;
        foreach(i; 0 .. massf.length) { massf[i] = w0*gs0.massf[i] + w1*gs1.massf[i]; }
        foreach(i; 0 .. rho_s.length) { rho_s[i] = w0*gs0.rho_s[i] + w1*gs1.rho_s[i]; }
        quality = w0*gs0.quality + w1*gs1.quality;
    }

    void copy_average_values_from(in GasState*[] others, GasModel gm)
    // Note that we must not send the current object in the others list as well.
    {
        size_t n = others.length;
        if (n == 0) throw new Error("Need to average from a nonempty array.");
        foreach(other; others) {
            if ( this is *other ) { throw new Error("Must not include destination in source list."); }
        }
        // Accumulate from a clean slate and then divide.
        p = 0.0; T = 0.0; u = 0.0; p_e = 0.0; a = 0.0;
        foreach(ref elem; u_modes) { elem = 0.0; }
        foreach(ref elem; T_modes) { elem = 0.0; }
        mu = 0.0; k = 0.0;
        foreach(ref elem; k_modes) { elem = 0.0; }
        sigma = 0.0;
        foreach(ref elem; massf) { elem = 0.0; }
        foreach(ref elem; rho_s) { elem = 0.0; }
        quality = 0.0;
        foreach(other; others) {
            p += other.p; T += other.T; u += other.u; p_e += other.p_e; a += other.a;
            foreach(i; 0 .. u_modes.length) { u_modes[i] += other.u_modes[i]; }
            foreach(i; 0 .. T_modes.length) { T_modes[i] += other.T_modes[i]; }
            mu += other.mu; k += other.k;
            foreach(i; 0 .. k_modes.length) { k_modes[i] += other.k_modes[i]; }
            sigma += other.sigma;
            foreach(i; 0 .. massf.length) { massf[i] += other.massf[i]; }
            foreach(i; 0 .. rho_s.length) { rho_s[i] += other.rho_s[i]; }
            quality += other.quality;
        }
        p /= n; T /= n; u /= n; p_e /= n; a /= n;
        foreach(ref elem; T_modes) { elem /= n; }
        foreach(ref elem; u_modes) { elem /= n; }
        mu /= n; k /= n;
        foreach(ref elem; k_modes) { elem /= n; }
        sigma /= n;
        foreach(ref elem; massf) { elem /= n; }
        foreach(ref elem; rho_s) { elem /= n; }
        quality /= n;
        // Now, make the properties consistent using the gas model.
        gm.update_thermo_from_pT(this);
        gm.update_sound_speed(this);
        gm.update_trans_coeffs(this);
    } // end copy_average_values_from()

    @nogc bool check_values(bool print_message=true) const
    {
        double RHOMIN = 0.0;
        double TMIN = 0.0;
        bool is_data_valid = true;
        if (!(isFinite(rho.re)) || rho < 1.01 * RHOMIN) {
            debug { if (print_message) writeln("Density invalid: ", rho); }
            is_data_valid = false;
        }
        if (!isFinite(T.re) || T < 1.01 * TMIN) {
            debug { if (print_message) writeln("Temperature invalid: ", T); }
            is_data_valid = false;
        }
        auto nmodes = u_modes.length;
        foreach(imode; 0 .. nmodes) {
            if (!isFinite(T_modes[imode].re) || T_modes[imode] < 1.01 * TMIN) {
                debug { if (print_message) writeln("Temperature[", imode, "] invalid: ", T_modes[imode]); }
                is_data_valid = false;
            }
            if ( !isFinite(u_modes[imode].re) ) {
                debug { if (print_message) writeln("Energy[", imode, "] invalid: ", u_modes[imode]); }
                is_data_valid = false;
            }
        }
        if (!isFinite(p.re)) {
            debug { if (print_message) writeln("Total pressure invalid: ", p); }
            is_data_valid = false;
        }
        if (!isFinite(p_e.re)) {
            debug { if (print_message) writeln("Electron pressure invalid: ", p_e); }
            is_data_valid = false;
        }
        if (!isFinite(a.re)) {
            debug { if (print_message) writeln("Sound speed invalid: ", a); }
            is_data_valid = false;
        }
        number f_sum = 0.0; foreach(elem; massf) f_sum += elem;
        if (f_sum < 0.99 || f_sum > 1.01 || !isFinite(f_sum.re)) {
            debug { if (print_message) writeln("Mass fraction sum bad: ", f_sum); }
            is_data_valid = false;
        }
        return is_data_valid;
    } // end check_values()

    string toString() const
    {
        char[] repr;
        repr ~= "GasState(";
        repr ~= "rho=" ~ to!string(rho);
        repr ~= ", p=" ~ to!string(p);
        repr ~= ", T=" ~ to!string(T);
        repr ~= ", u=" ~ to!string(u);
        repr ~= ", p_e=" ~ to!string(p_e);
        repr ~= ", a=" ~ to!string(a);
        repr ~= ", T_modes=" ~ to!string(T_modes);
        repr ~= ", u_modes=" ~ to!string(u_modes);
        repr ~= ", mu=" ~ to!string(mu);
        repr ~= ", k=" ~ to!string(k);
        repr ~= ", k_modes=" ~ to!string(k_modes);
        repr ~= ", massf=" ~ to!string(massf);
        repr ~= ", rho_s=" ~ to!string(rho_s);
        repr ~= ", quality=" ~ to!string(quality);
        repr ~= ", sigma=" ~ to!string(sigma);
        repr ~= ")";
        return to!string(repr);
    }

version(complex_numbers) {
    @nogc
    void clear_imaginary_components()
    // When performing the complex-step Frechet derivative in the Newton-Krylov accelerator,
    // the flowstate values accumulate imaginary components, so we have to start with a clean slate, so to speak.
    {
        rho.im = 0.0;
        p.im = 0.0;
        p_e.im = 0.0;
        a.im = 0.0;
        T.im = 0.0;
        u.im = 0.0;
        foreach( i; 0..u_modes.length) u_modes[i].im = 0.0;
        foreach( i; 0..T_modes.length) T_modes[i].im = 0.0;
        mu.im = 0.0;
        k.im = 0.0;
        foreach( i; 0..k_modes.length) k_modes[i].im = 0.0;
        sigma.im = 0.0;
        foreach( i; 0..massf.length) massf[i].im = 0.0;
        foreach( i; 0..rho_s.length) rho_s[i].im = 0.0;
        quality.im = 0.0;
    } // end clear_imaginary_components()
} // end version(complex)

} // end class GasState
