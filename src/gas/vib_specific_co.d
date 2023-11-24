/*
    Prototyping for vibrationally specific CO gas model, intended for use in simulating 
    gas dynamic lasers.

    Notes:
    "Kinetic Modelling of the High-Power Carbon Monoxide Laser*"
    Joseph W. Rich, Journal of Applied Physics, Volume 42, Number 7, June 1971

    @author: Nick Gibbons (24/03/22)
*/

module gas.vib_specific_co;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;

import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.composite_gas;

immutable uint _n_vibe_states = 10;

class VibSpecificCO: GasModel {
    this(lua_State *L)
    {
        type_str = "VibSpecificCO";

        // Set up the species names array
        _n_modes = 0;
        _n_species = _n_vibe_states;
        _species_names.length = _n_species;
        foreach(i; 0 .. _n_vibe_states){
            _species_names[i] = format("CO-%02d", i);
        }
        create_species_reverse_lookup();

        // Set up molar masses, including possible extra species
        _mol_masses.length = _n_species;
        foreach (isp; 0 .. _n_vibe_states) {
            _mol_masses[isp] = _M;
        }

        // The energy levels of the individual vibration levels are constant, once computed.
        _vib_energy.length = _n_vibe_states;
        _vib_energy_per_kg.length = _n_vibe_states;
        foreach (i; 0 .. _n_vibe_states) {
            _vib_energy[i] = vib_energy(i);
            _vib_energy_per_kg[i]= (Avogadro_number/_M) * _vib_energy[i];
        }

        // TODO: Viscous stuff
        //gs = GasState(1, 0); // Fake gas state for interfacing with transport properties.
        //mTransProps = new GasMixtureTransProps(L, ["CO"]);
    }

    this(string gas_file_name){
        lua_State* L = init_lua_State();
        doLuaFile(L, gas_file_name);
        this(L);
        lua_close(L);
    }

    override string toString() const
    {
        char[] repr;
        repr ~= format("VibSpecificCO=(numVibLevels=%d)", _n_vibe_states);
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        Q.rho = Q.p/(_R*Q.T);
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_rhou(ref GasState Q)
    {
        // From internal energy, remove vibrational energy before computing trans-rotational temperature.
        number u = Q.u;
        foreach (i; 0 .. _n_vibe_states) { u -= Q.massf[i]*_vib_energy_per_kg[i]; }
        Q.T = (0.4/_R)*u;
        Q.p = Q.rho*_R*Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q)
    {
        Q.p = Q.rho*_R*Q.T;
        // Start with trans-rotational component of internal energy and add vibrational energy.
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_rhop(ref GasState Q)
    {
        Q.T = Q.p/(Q.rho*_R);
        // Start with trans-rotational component of internal energy and add vibrational energy.
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        number sum=0.0;
        foreach(v; 0 .. _n_vibe_states) sum += Q.massf[v]*log(Q.massf[v]);
        Q.T = Tref*exp((s+_R*log(Q.p/pref) + _R*sum)/_cp);
        Q.rho = Q.p/(_R*Q.T);
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        throw new Error("VibSpecificCO:update_thermo_from_hs NOT IMPLEMENTED.");
    }
    override void update_sound_speed(ref GasState Q)
    {
        Q.a = sqrt(_gamma*_R*Q.T);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        // FIXME: Put the fake gas state stuff here
        //mTransProps.updateTransProps(gs);
        //throw new Error("ERROR: VibSpecificCO:update_trans_coeffs NOT IMPLEMENTED.");
        Q.mu = 0.0; Q.k = 0.0;
    }
    override number dudT_const_v(in GasState Q)
    {
        return to!number(_cv);
    }
    override number dhdT_const_p(in GasState Q)
    {
        return to!number(_cp);
    }
    override number dpdrho_const_T(in GasState Q)
    {
        return _R*Q.T;
    }
    override number gas_constant(in GasState Q)
    {
        return to!number(_R);
    }
    override number internal_energy(in GasState Q)
    {
        number u = 0.0;
        number cvT = _cv*Q.T;
        foreach (i; 0 .. _n_vibe_states) { u += Q.massf[i]*(cvT + _vib_energy_per_kg[i]); }
        return u;
    }
    override number enthalpy(in GasState Q)
    {
        number h = 0.0;
        number cpT = _cp*Q.T;
        foreach (i; 0 .. _n_vibe_states) { h += Q.massf[i]*(cpT + _vib_energy_per_kg[i]); }
        return h;
    }
    override number entropy(in GasState Q)
    {
    /*
        Not completely sure about this. See 06/04/22 derivation. I think I can get away with
        setting s @ 300K and 1 atm to zero, which eliminates the reference entropy.
    */
        number s=0.0;
        foreach(i; 0 .. n_vibe_states){
            number Xi = Q.massf[i]; // Molefraction = mass fraction in CO only mode
            number si = _cp*log(Q.T/Tref) - _R*log(Xi*Q.p/pref);
            s += Q.massf[i]*si;
        }
        return s;
    }

    @nogc uint n_vibe_states() const { return _n_vibe_states; }

    @nogc const
    number boltzmann_eq_population_fraction(size_t v, number T){
    /*
        Returns the equilibrium population fraction for quantum-level v, given temperature.
        v==0 is the ground state.
    */
        if (v >= _n_vibe_states) return to!number(0.0);

        number summ = 0.0;
        foreach(ej; 0 .. _n_vibe_states) { summ += exp(-_vib_energy[ej] / (Boltzmann_constant*T)); }
        number nvf = exp(-_vib_energy[v]/(Boltzmann_constant*T)) / summ;
        return nvf;
    }
protected:
    // General Gas Stuff
    //GasState gs;
    //TransportPropertiesModel mTransProps;
    immutable double Tref = 300.0;
    immutable double pref = 101.325e3;
    immutable double Ru = 8.31446261815324; // Exact as of 2019 SI redefinition
    immutable double _M = 0.0280101;
    immutable double _R = Ru/_M;
    immutable double _Cv = 5.0/2.0*Ru;
    immutable double _Cp = _Cv + 1.0*Ru;
    immutable double _cv = _Cv/_M;
    immutable double _cp = _Cp/_M;
    immutable double _gamma = _cp/_cv;

    // State Specific Constants
    immutable double E01 = 4.31e-13*1e-7; // erg -> J Ground state energy (FIXME: Is this E1? A typo???)
    immutable double d   = 0.00598;       //          CO Anharmonicity

    // Internal storage
    double[] _vib_energy; // quantum level energies
    double[] _vib_energy_per_kg;

    @nogc const double vib_energy(uint v){
    // Returns the quantum level energy for mode v, subtracting off the zero-point so that E(v=0)=0
        return E01*(v - d*v*(v+1.0));
    }

    // There used to be a Tvib function here but I don't think we need it
} // end class VibSpecificCO

class VibSpecificCOMixture: VibSpecificCO {
    // Public class specific properties that the kinetics model needs
    uint n_others;
    GasState cgs;
    CompositeGas cgm;

    this(lua_State *L)
    {
        super(L);
        type_str = "VibSpecificCOMixture";

        // We have additional species, which are not treated state specifically
        string otherSpeciesFile = getString(L, "other_species");
        cgm = new CompositeGas(otherSpeciesFile);
        n_others = cgm.n_species;
        cgs = GasState(cgm);

        // Set up the species names array
        _n_species = _n_vibe_states + n_others;
        _species_names.length = _n_species;
        foreach(i; 0 .. n_others){
            _species_names[_n_vibe_states+i] = cgm.species_name(i);
        }
        create_species_reverse_lookup();

        // Add the extra species molar masses to the end of the _mol_masses array
        _mol_masses.length = _n_species;
        foreach(isp; 0 .. n_others){
            _mol_masses[_n_vibe_states+isp] = cgm.mol_masses[isp];
        }
    }

    this(string gas_file_name){
        lua_State* L = init_lua_State();
        doLuaFile(L, gas_file_name);
        this(L);
        lua_close(L);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        number R = gas_constant(Q);
        Q.rho = Q.p/R/Q.T;
        Q.u = internal_energy(Q);
    }

    override void update_thermo_from_rhou(ref GasState Q)
    {
        Q.T = temperature_from_energy(Q);
        number R = gas_constant(Q);
        Q.p = Q.rho*R*Q.T;
    }

    override void update_thermo_from_rhoT(ref GasState Q)
    {
        number R = gas_constant(Q);
        Q.p = Q.rho*R*Q.T;
        Q.u = internal_energy(Q);
    }

    override void update_thermo_from_rhop(ref GasState Q)
    {
        number R = gas_constant(Q);
        Q.T = Q.p/R/Q.rho;
        Q.u = internal_energy(Q);
    }

    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        throw new Error("VibSpecificCOMixture:update_thermo_from_ps NOT IMPLEMENTED.");
    }

    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        throw new Error("VibSpecificCOMixture:update_thermo_from_hs NOT IMPLEMENTED.");
    }

    override void update_sound_speed(ref GasState Q)
    {
        number R = gas_constant(Q);
        number gamma = dhdT_const_p(Q)/dudT_const_v(Q);
        Q.a = sqrt(gamma*R*Q.T);
    }

    override void update_trans_coeffs(ref GasState Q)
    {
        // FIXME: This is going to be really nasty
        Q.mu = 0.0; Q.k = 0.0;
    }

    override number dudT_const_v(in GasState Q)
    {
        number YCO = massf_of_CO(Q);
        set_cgs_massf(Q);
        cgs.T = Q.T;
        return YCO*_cv + cgm.dudT_const_v(cgs);
    }

    override number dhdT_const_p(in GasState Q)
    {
        number YCO = massf_of_CO(Q);
        set_cgs_massf(Q);
        cgs.T = Q.T;
        return YCO*_cp + cgm.dhdT_const_p(cgs);
    }

    override number dpdrho_const_T(in GasState Q)
    {
        number R = gas_constant(Q);
        return R*Q.T;
    }

    override number gas_constant(in GasState Q)
    {
        number YCO = massf_of_CO(Q);
        set_cgs_massf(Q);
        number Ro = cgm.gas_constant(cgs);
        return YCO*_R + Ro;
    }

    override number internal_energy(in GasState Q)
    {
        number uCO = super.internal_energy(Q);
        set_cgs_massf(Q);
        cgs.T = Q.T;
        number uo = cgm.internal_energy(cgs);
        return uCO+uo;
    }

    override number enthalpy(in GasState Q)
    {
        number hCO = super.enthalpy(Q);
        set_cgs_massf(Q);
        cgs.T = Q.T;
        number ho = cgm.enthalpy(cgs);
        return hCO+ho;
    }

    override number entropy(in GasState Q)
    {
        throw new Error("VibSpecificCOMixture:entropy NOT IMPLEMENTED.");
    }


private:
    @nogc
    number temperature_from_energy(ref GasState gs)
    {
        /*
            Newton's method that is very similar to the one in gas/thermo/two_temperature_gas.d
        */
        int MAX_ITERATIONS = 20;

        // Take the supplied T as the initial guess.
        number T_guess = gs.T;
        number u0 = internal_energy(gs);
        number f_guess =  u0 - gs.u;

        // Before iterating, check if the supplied guess is good enough.
        double E_TOL = 1e-5;
        if (fabs(f_guess) < E_TOL) {

            version(complex_numbers) {
            /*
                Set complex components using the analytical derivatives
            */
                number Cv = dudT_const_v(gs); // I don't think this can go bad
                set_cgs_massf(gs);
                cgs.T = T_guess;

                // The analytical derivatives involve Cv and the species energies
                double T_guess_im = gs.u.im;
                foreach(i; 0 .. _n_vibe_states) T_guess_im -= gs.massf[i].im*(_cv*T_guess.re + _vib_energy_per_kg[i]);
                foreach(i; 0 .. n_others)       T_guess_im -= gs.massf[i+_n_vibe_states].im*cgm.energyPerSpeciesInMode(cgs, i, 0).re;
                T_guess_im /= Cv.re;
                T_guess.im = T_guess_im;
            }
            return T_guess;
        }

        // We'll keep adjusting our temperature estimate
        // until it is less than TOL.
        double TOL = 1.0e-9;

        // Begin iterating.
        int count = 0;
        number Cv, dT;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            Cv = dudT_const_v(gs);
            dT = -f_guess/Cv;
            T_guess += dT;
            gs.T = T_guess;
            if (fabs(dT) < TOL) {
                break;
            }
            f_guess = internal_energy(gs) - gs.u;
            count++;
        }

        if ((count == MAX_ITERATIONS)&&(fabs(dT)>1e-3)) {
            string msg = "The 'temperature_from_energy' function failed to converge.\n";
            debug {
                msg ~= format("The final value for T was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= gs.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }
        return T_guess;
    }

    @nogc
    number massf_of_CO(in GasState gs) const {
        number YCO = 0.0;
        foreach(i; 0 .. _n_vibe_states){
            YCO += gs.massf[i];
        }
        return YCO;
    }

    @nogc
    void set_cgs_massf(in GasState gs) {
        foreach(i; 0 .. n_others){
            cgs.massf[i] = gs.massf[_n_vibe_states+i];
        }
    }
}

version(vib_specific_co_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        //lua_State* L = init_lua_State();
        //doLuaFile(L, );
        //auto gm = new VibSpecificCO(L);
        //lua_close(L);
        auto gm = new VibSpecificCO("sample-data/vib-specific-CO-gas.lua");
        auto Q = GasState(gm.n_species, 0);

        // Practice problem to match the data in table II
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        // Set up the species mass fractions assuming equilibrium.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = gm.boltzmann_eq_population_fraction(v, Q.T); }

        double R_CO = 296.8379191791532; // gas constant for CO
        double M_CO = 0.0280101; // kg/mole
        double gamma = 7.0/5.0; // ratio of specific heats.

        gm.update_thermo_from_pT(Q);
        double my_rho = 26.7 / (R_CO * 175.0);
        assert(isClose(Q.rho, my_rho, 1.0e-6), failedUnitTest());

        double my_u = 2.5 * R_CO * 175.0;
        foreach (i; 0 .. gm.n_vibe_states) {
            my_u += (Avogadro_number/M_CO) * gm.vib_energy(i) * Q.massf[i];
        }
        assert(isClose(Q.u, my_u, 1.0e-6), failedUnitTest());

        //gm.update_trans_coeffs(Q);
        //assert(isClose(Q.mu, 0.0, 1.0e-6), failedUnitTest());
        //assert(isClose(Q.k, 0.0, 1.0e-6), failedUnitTest());

        gm.update_sound_speed(Q);
        double my_a = sqrt(gamma * R_CO * 175.0);
        assert(isClose(Q.a, my_a, 1.0e-6), failedUnitTest());

        return 0;
    }
}

version(vib_specific_co_mixture_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        auto gm = new VibSpecificCOMixture("sample-data/vib-specific-CO-mixture.lua");
        auto Q = GasState(gm.n_species, 0);

        // Practice problem to match the data in table II
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        // Set up the species mass fractions assuming equilibrium and no N2.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = gm.boltzmann_eq_population_fraction(v, Q.T); }
        Q.massf[gm.n_vibe_states] = 0.0;

        double R_CO = 296.8379191791532; // gas constant for CO
        double M_CO = 0.0280101; // kg/mole
        double gamma = 7.0/5.0; // ratio of specific heats.

        gm.update_thermo_from_pT(Q);
        double my_rho = 26.7 / (R_CO * 175.0);
        assert(isClose(Q.rho, my_rho, 1.0e-6), failedUnitTest());

        double my_u = 2.5 * R_CO * 175.0;
        foreach (i; 0 .. gm.n_vibe_states) {
            my_u += (Avogadro_number/M_CO) * gm.vib_energy(i) * Q.massf[i];
        }
        assert(isClose(Q.u, my_u, 1.0e-6), failedUnitTest());

        //gm.update_trans_coeffs(Q);
        //assert(isClose(Q.mu, 0.0, 1.0e-6), failedUnitTest());
        //assert(isClose(Q.k, 0.0, 1.0e-6), failedUnitTest());

        gm.update_sound_speed(Q);
        double my_a = sqrt(gamma * R_CO * 175.0);
        assert(isClose(Q.a, my_a, 1.0e-6), failedUnitTest());

        number save_T = Q.T;
        Q.T = 300.0;
        gm.update_thermo_from_rhou(Q);
        assert(isClose(Q.T, save_T, 1.0e-6), failedUnitTest());

        // ---------------------------------------------------------------------------
        // Now check that we match the N2 model when only N2 is present
        auto cgm = new CompositeGas("sample-data/thermally-perfect-n2.lua");
        auto cgs = GasState(cgm);
        cgs.p = 26.7;  // Pa
        cgs.T = 175.0; // K
        cgs.massf[0] = 1.0;
        cgm.update_thermo_from_pT(cgs);
        cgm.update_sound_speed(cgs);

        // Practice problem to match the data in table II
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        // Set up the species mass fractions assuming equilibrium and no N2.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = 0.0; }
        Q.massf[gm.n_vibe_states] = 1.0;
        gm.update_thermo_from_pT(Q);
        gm.update_sound_speed(Q);

        assert(isClose(Q.rho, cgs.rho, 1.0e-6), failedUnitTest());
        assert(isClose(Q.u, cgs.u, 1.0e-6), failedUnitTest());
        assert(isClose(Q.a, cgs.a, 1.0e-6), failedUnitTest());

        save_T = Q.T;
        Q.T = 300.0;
        gm.update_thermo_from_rhou(Q);
        assert(isClose(Q.T, save_T, 1.0e-6), failedUnitTest());

        return 0;
    }
}
