/**
 * therm_perf_gas_mix.d
 * Thermally perfect gas mix model for use in the CFD codes.
 *
 * Author: Rowan G. and Peter J.
 * First code: 27-Jan-2015
 */

module gas.therm_perf_gas;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import std.algorithm;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;
import nm.brent;
import nm.bracketing;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.thermo.cea_thermo_curves;
import gas.thermo.perf_gas_mix_eos;
import gas.thermo.therm_perf_gas_mix_eos;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.cea_viscosity;
import gas.diffusion.cea_therm_cond;
import gas.diffusion.chemkin_viscosity;
import gas.diffusion.chemkin_therm_cond;
import gas.diffusion.sutherland_viscosity;
import gas.diffusion.sutherland_therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;
import gas.diffusion.gasgiant_transport_properties;

class ThermallyPerfectGas: GasModel {
public:
    this(lua_State* L)
    // Construct the model from parameters that are contained in a Lua interpreter.
    {
        type_str = "ThermallyPerfectGas";
        getArrayOfStrings(L, "species", _species_names);
        _n_species = cast(uint) _species_names.length;
        _n_modes = 0; // Single temperature gas
        if (canFind(_species_names, "e-") || canFind(_species_names, "eminus")) {
            if (!(_species_names[$-1] == "e-" || _species_names[$-1] == "eminus")) {
                throw new Error("Electrons should be last species.");
            }
            _is_plasma = true;
        }
        bool useGasGiantTransProps = false;
        lua_getglobal(L, "use_gas_giant_transport_properties");
        if (!lua_isnil(L, -1)) useGasGiantTransProps = true;
        lua_pop(L, 1);

        if (useGasGiantTransProps) {
            // Check early we have correct number of species for using gas giant model
            if (_n_species != 6) {
                string errMsg = "Incorrect number of species to use the gas giant transport properties model.\n";
                errMsg ~= format("6 species are expected, but %d species have been given.\n", _n_species);
                throw new Error(errMsg);
            }
        }

        // 0. Initialise private work arrays
        _Cp.length = _n_species;
        _Cv.length = _n_species;
        _h.length = _n_species;
        _s.length = _n_species;
        _molef.length = _n_species;
        // 1. Initialise gas constants from molecular mass
        //    and grab charge.
        _R.length = _n_species;
        _mol_masses.length = _n_species;
        _charge.length = n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _mol_masses[isp] = getDouble(L, -1, "M");
            _R[isp] = R_universal/_mol_masses[isp];
            _charge[isp] = getInt(L, -1, "charge");
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        create_species_reverse_lookup();
        // 1b. Gather L-J parameters
        _LJ_sigmas.length = _n_species;
        _LJ_epsilons.length = _n_species;
        _Lewis_numbers.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _LJ_sigmas[isp] = getDouble(L, -1, "sigma");
            _LJ_epsilons[isp] = getDouble(L, -1, "epsilon");
            _Lewis_numbers[isp] = getDouble(L, -1, "Lewis");
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        // 2. Set the p-v-T EOS
        _pgMixEOS = new PerfectGasMixEOS(_R);
        // 3. Set the e-v-T EOS (which in this case is an e-T relationship only)
        foreach ( isp; 0.._n_species ) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            lua_getfield(L, -1, "thermoCoeffs");
            _curves ~= new CEAThermoCurve(L, _R[isp]);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        _tpgMixEOS = new ThermallyPerfectGasMixEOS(_R, _curves);
        // 4. Set the viscosity model
        if (useGasGiantTransProps) {
            // Before we go on, we need to check the species are in correct order.
            bool speciesOrderCorrect = true;
            if (species_index("H2") != gas.diffusion.gasgiant_transport_properties.Sp.H2) speciesOrderCorrect = false;
            if (species_index("H") != gas.diffusion.gasgiant_transport_properties.Sp.H) speciesOrderCorrect = false;
            if (species_index("H+") != gas.diffusion.gasgiant_transport_properties.Sp.Hp) speciesOrderCorrect = false;
            if (species_index("He") != gas.diffusion.gasgiant_transport_properties.Sp.He) speciesOrderCorrect = false;
            if (species_index("He+") != gas.diffusion.gasgiant_transport_properties.Sp.Hep) speciesOrderCorrect = false;
            if ((species_index("e-") != gas.diffusion.gasgiant_transport_properties.Sp.e) &&
                (species_index("eminus") != gas.diffusion.gasgiant_transport_properties.Sp.e))
                speciesOrderCorrect = false;
            if (!speciesOrderCorrect) {
                string errMsg = "The species have been declared in the wrong order if using the gas giant transport properties model.";
                errMsg ~= "The correct order is: H2, H, H+, He, He+, e-\n";
                throw new Error(errMsg);
            }
            string dbOption = getString(L, "gas_giant_visc_db");
            _viscModel = new GasGiantViscosity(dbOption);
        }
        else {

            Viscosity[] vms;
            foreach ( isp; 0.._n_species ) {
                lua_getglobal(L, "db");
                lua_getfield(L, -1, _species_names[isp].toStringz);
                lua_getfield(L, -1, "viscosity");
                lua_getfield(L, -1, "model");
                string model = to!string(luaL_checkstring(L, -1));
                lua_pop(L, 1);
                switch (model) {
                case "CEA":
                    vms ~= createCEAViscosity(L);
                    break;
                case "Sutherland":
                    vms ~= createSutherlandViscosity(L);
                    break;
                case "Chemkin":
                    vms ~= createChemkinViscosity(L);
                    break;
                default:
                    string errMsg = format("The viscosity model '%s' is not available.\n", model);
                    errMsg ~= format("This error occurred for species %d when constructing "
                                     ~"the thermally perfect gas mix.\n", isp);
                    throw new Error(errMsg);
                }
                lua_pop(L, 1);
                lua_pop(L, 1);
                lua_pop(L, 1);
            }
            _viscModel = new WilkeMixingViscosity(vms, _mol_masses);
        }
        // 5. Set the thermal conductivity model
        if (useGasGiantTransProps) {
            string dbOption = getString(L, "gas_giant_tc_db");
            string modelOption = getString(L, "gas_giant_tc_model");
            _thermCondModel = new  GasGiantThermalConductivity(dbOption, modelOption);
        }
        else {
            ThermalConductivity[] tcms;
            foreach (isp; 0.._n_species ) {
                lua_getglobal(L, "db");
                lua_getfield(L, -1, _species_names[isp].toStringz);
                lua_getfield(L, -1, "thermal_conductivity");
                lua_getfield(L, -1, "model");
                string model = to!string(luaL_checkstring(L, -1));
                lua_pop(L, 1);
                switch (model) {
                case "CEA":
                    tcms ~= createCEAThermalConductivity(L);
                    break;
                case "Sutherland":
                    tcms ~=  createSutherlandThermalConductivity(L);
                    break;
                case "Chemkin":
                    tcms ~=  createChemkinThermalConductivity(L);
                    break;
                default:
                    string errMsg = format("The thermal conductivity model '%s' "
                                           ~"is not available.\n", model);
                    errMsg ~= format("This error occurred for species %d when constructing "
                                     ~"the thermally perfect gas mix.\n", isp);
                    throw new Error(errMsg);
                }
                lua_pop(L, 1);
                lua_pop(L, 1);
                lua_pop(L, 1);
            }
            _thermCondModel = new WilkeMixingThermCond(tcms, _mol_masses);
        }
    } // end constructor using Lua interpreter

    void attachGasModelToGasGiantModel()
    {
        GasGiantViscosity ggv = cast(GasGiantViscosity) _viscModel;
        ggv.attachGasModel(this);
        GasGiantThermalConductivity ggtc = cast(GasGiantThermalConductivity) _thermCondModel;
        ggtc.attachGasModel(this);
    }

    this(in string fname)
    {
        auto L = init_lua_State();
        doLuaFile(L, fname);
        this(L);
        lua_close(L); // We no longer need the Lua interpreter.
    } // end constructor from a Lua file

    override string toString() const
    {
        return "ThermallyPerfectGas(species=[TODO])";
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        _pgMixEOS.update_density(Q);
        _tpgMixEOS.update_energy(Q);
    }
    override void update_thermo_from_rhou(ref GasState Q)
    {
        _tpgMixEOS.update_temperature(Q);
        _pgMixEOS.update_pressure(Q);
    }
    override void update_thermo_from_rhoT(ref GasState Q)
    {
        _tpgMixEOS.update_energy(Q);
        _pgMixEOS.update_pressure(Q);
    }
    override void update_thermo_from_rhop(ref GasState Q)
    {
        _pgMixEOS.update_temperature(Q);
        _tpgMixEOS.update_energy(Q);
    }
    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        double TOL = 1.0e-6;
        number delT = 100.0;
        number T1 = Q.T;
        // It's possible that T is 'nan' if it's never been set, so:
        if ( isNaN(T1) )
            T1 = 300.0; // Just set at some value to get started.
        number Tsave = T1;
        number T2 = T1 + delT;

        number zeroFun(number T) {
            Q.T = T;
            number s_guess = entropy(Q);
            return s - s_guess;
        }

        number T_lower_limit = 10.0;
        if ( bracket!(zeroFun,number)(T1, T2, T_lower_limit) == -1 ) {
            string msg = "The 'bracket' function failed to find temperature values\n";
            debug {
                msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_ps().\n";
                msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
            }
            throw new Exception(msg);
        }

        if ( T1 < T_MIN )
            T1 = T_MIN;

        try {
            Q.T = solve!(zeroFun,number)(T1, T2, TOL);
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find temperature\n";
            debug {
                msg ~= "in function ThermallyPerfectGas.update_thermo_from_ps().\n";
                msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                msg ~= format("The target entropy value was: %12.6f\n", s);
                msg ~= format("The GasState is currently:\n");
                msg ~= Q.toString();
                msg ~= "The message from the ridder.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }
        _tpgMixEOS.update_energy(Q);
        _pgMixEOS.update_density(Q);
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        // We do this in two stages.
        // First, from enthalpy we compute temperature.
        double TOL = 1.0e-6;
        number delT = 100.0;
        number T1 = Q.T;
        // It's possible that T is 'nan' if it's never been set, so:
        if ( isNaN(T1) )
            T1 = 300.0; // Just set at some value to get started.
        number Tsave = T1;
        number T2 = T1 + delT;

        number zeroFun(number T)
        {
            Q.T = T;
            number h_guess = enthalpy(Q);
            return h - h_guess;
        }

        if ( bracket!(zeroFun,number)(T1, T2) == -1 ) {
            string msg = "The 'bracket' function failed to find temperature values";
            debug {
                msg ~= "\nthat bracketed the zero function in ThermallyPerfectGas.update_thermo_from_hs().\n";
                msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
            }
            throw new Exception(msg);
        }

        if ( T1 < T_MIN )
            T1 = T_MIN;

        try {
            Q.T = solve!(zeroFun,number)(T1, T2, TOL);
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find temperature";
            debug {
                msg ~= "\nin function ThermallyPerfectGas.update_thermo_from_hs().\n";
                msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                msg ~= format("The target enthalpy value was: %12.6f\n", h);
                msg ~= format("The GasState is currently:\n");
                msg ~= Q.toString();
                msg ~= "The message from the ridder.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }

        // Second, we can iterate to find the pressure that gives
        // correct entropy.
        TOL = 1.0e-3;
        number delp = 1000.0;
        number p1 = Q.p;
        number psave = p1;
        number p2 = p1 + delp;

        number zeroFun2(number p)
        {
            Q.p = p;
            number s_guess = entropy(Q);
            return s - s_guess;
        }

        if ( bracket!(zeroFun2,number)(p1, p2) == -1 ) {
            string msg = "The 'bracket' function failed to find pressure values";
            debug {
                msg ~= "\nthat bracketed the zero function in ThermallyPerfectGas.update_thermo_from_hs().\n";
                msg ~= format("The final values are: p1 = %12.6f and p2 = %12.6f\n", p1, p2);
            }
            throw new Exception(msg);
        }

        if ( p1 < 0.0 )
            p1 = 1.0;

        try {
            Q.p = solve!(zeroFun2,number)(p1, p2, TOL);
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find pressure";
            debug {
                msg ~= "\nin function ThermallyPerfectGas.update_thermo_from_hs().\n";
                msg ~= format("The initial pressure guess was: %12.6f\n", psave);
                msg ~= format("The target entropy value was: %12.6f\n", s);
                msg ~= format("The GasState is currently:\n");
                msg ~= Q.toString();
                msg ~= "The message from the ridder.solve function is:\n";
                msg ~= e.msg;
            }
            throw new Exception(msg);
        }
        _tpgMixEOS.update_energy(Q);
        _pgMixEOS.update_density(Q);
    }
    override void update_sound_speed(ref GasState Q)
    {
        // Reference:
        // Cengel and Boles (1998)
        // Thermodynamics: an Engineering Approach, 3rd edition
        // McGraw Hill
        // Equation 16-10 on p. 849

        // "frozen" sound speed
        Q.a = sqrt(gamma(Q)*dpdrho_const_T(Q));
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        _viscModel.update_viscosity(Q);
        _thermCondModel.update_thermal_conductivity(Q);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
        throw new Exception("not implemented");
    }
    */
    override number dudT_const_v(in GasState Q)
    {
        // Noting that Cv = Cp - R
        foreach ( i; 0.._n_species ) { _Cv[i] = _curves[i].eval_Cp(Q.T) - _R[i]; }
        return mass_average(Q, _Cv);
    }
    override number dhdT_const_p(in GasState Q)
    {
        foreach ( i; 0.._n_species ) { _Cp[i] = _curves[i].eval_Cp(Q.T); }
        return mass_average(Q, _Cp);
    }
    override number dpdrho_const_T(in GasState Q)
    {
        number R = gas_constant(Q);
        return R*Q.T;
    }
    override number gas_constant(in GasState Q)
    {
        return mass_average(Q, _R);
    }
    override number gas_constant(in GasState Q, int isp)
    {
        return to!number(_R[isp]);
    }
    override number internal_energy(in GasState Q)
    {
        return Q.u;
    }
    override number internal_energy(in GasState Q, int isp)
    {
        number h = _curves[isp].eval_h(Q.T);
        number u = h - _R[isp]*Q.T;
        return u;
    }
    override number enthalpy(in GasState Q)
    {
        foreach ( isp; 0.._n_species) {
            _h[isp] = _curves[isp].eval_h(Q.T);
        }
        return mass_average(Q, _h);
    }
    override number enthalpy(in GasState Q, int isp)
    {
        return _curves[isp].eval_h(Q.T);
    }
    override number entropy(in GasState Q)
    {
        foreach ( isp; 0.._n_species ) {
            _s[isp] = _curves[isp].eval_s(Q.T) - _R[isp]*log(Q.p/P_atm);
        }
        return mass_average(Q, _s);
    }
    override number entropy(in GasState Q, int isp)
    {
        return _curves[isp].eval_s(Q.T);
    }

    override number Cp(in GasState Q, int isp)
    {
        return _curves[isp].eval_Cp(Q.T);
    }

    override void balance_charge(ref GasState Q)
    {
        if (_is_plasma) {
            massf2molef(Q, _molef);
            number molefIons = 0.0;
            // Loop to n_species - 1 so that we do NOT include electron
            // in counting up ionic contributions.
            foreach (isp; 0 .. _n_species-1) {
                molefIons += _charge[isp] * _molef[isp];
            }
            // Now set electron mole fraction to balance
            _molef[_n_species-1] = molefIons;
            molef2massf(_molef, Q);
        }
    }
    @nogc
    override void binary_diffusion_coefficients(ref const(GasState) Q, ref number[][] D)
    {
        // Expression from:
        // Reid et al.
        // The Properties of Gases and Liquids

        // Loop through only the upper elements.
        // Moved here by NNG 09/09/2020
        debug{ assert(D.length==_n_species); }
        number T = Q.T;
        number p = Q.p;
        foreach (isp; 0 .. _n_species) {
            foreach (jsp; isp+1 .. _n_species) {
                // These three lines can be precomputed to save some time.
                number _eps = sqrt(_LJ_epsilons[isp] * _LJ_epsilons[jsp]);
                number _M = 2.0/(1.0/_mol_masses[isp] + 1.0/_mol_masses[jsp])*1.0e3;
                number _sigma = 0.5*(_LJ_sigmas[isp] + _LJ_sigmas[jsp]);

                number T_star = T/_eps;
                number omega = _A/(pow(T_star, _B));
                omega += _C/(exp(_D*T_star));
                omega += _E/(exp(_F*T_star));
                omega += _G/(exp(_H*T_star));

                number denom = p/P_atm;
                denom *= sqrt(_M);
                denom *= _sigma*_sigma;
                denom *= omega;
                number numer = 0.00266*sqrt(T*T*T);
                numer *= 1.0e-4; // cm^2/s --> m^2/s
                D[isp][jsp] = numer/denom;
                D[jsp][isp] = D[isp][jsp];
            }
        }
    }


protected:
    double[] _R;
    PerfectGasMixEOS _pgMixEOS;
    ThermallyPerfectGasMixEOS _tpgMixEOS;
    CEAThermoCurve[] _curves;
    ViscosityMixtureModel _viscModel;
    ThermalConductivityMixtureModel _thermCondModel;
    // Working array space
    number[] _Cp, _Cv, _h, _s, _molef;
    // Data for binary diffusion calculations
    immutable double _A = 1.06036;
    immutable double _B = 0.15610;
    immutable double _C = 0.19300;
    immutable double _D = 0.47635;
    immutable double _E = 1.03587;
    immutable double _F = 1.52996;
    immutable double _G = 1.76474;
    immutable double _H = 3.89411;
} // end class ThermallyPerfectGas


version(therm_perf_gas_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
        //
        auto gm = new ThermallyPerfectGas("sample-data/therm-perf-5-species-air.lua");
        auto gd = GasState(5, 0);
        assert(isClose(3.621, gm.LJ_sigmas[0]), failedUnitTest());
        assert(isClose(97.530, gm.LJ_epsilons[0]), failedUnitTest());

        gd.p = 1.0e6;
        gd.T = 2000.0;
        gd.massf = [to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2)];
        gm.update_thermo_from_pT(gd);
        assert(approxEqualNumbers(to!number(11801825.6), gd.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.2840117), gd.rho, 1.0e-6), failedUnitTest());

        gd.rho = 2.0;
        gd.u = 14.0e6;
        gm.update_thermo_from_rhou(gd);
        assert(approxEqualNumbers(to!number(3373757.4), gd.p, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(4331.944), gd.T, 1.0e-6), failedUnitTest());

        gd.T = 10000.0;
        gd.rho = 1.5;
        gm.update_thermo_from_rhoT(gd);
        assert(approxEqualNumbers(to!number(5841068.3), gd.p, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(20340105.9), gd.u, 1.0e-6), failedUnitTest());

        gd.rho = 10.0;
        gd.p = 5.0e6;
        gm.update_thermo_from_rhop(gd);
        assert(approxEqualNumbers(to!number(11164648.5), gd.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1284.012), gd.T, 1.0e-6), failedUnitTest());

        gd.p = 1.0e6;
        number s = 10000.0;
        gm.update_thermo_from_ps(gd, s);
        assert(approxEqualNumbers(to!number(2560.118), gd.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(12313952.52), gd.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.00309), gd.rho, 1.0e-6), failedUnitTest());

        s = 11000.0;
        number h = 17.0e6;
        gm.update_thermo_from_hs(gd, h, s);
        assert(approxEqualNumbers(to!number(5273.103), gd.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(14946629.7), gd.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(0.4603513), gd.rho, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(945271.84), gd.p, 1.0e-4), failedUnitTest());

        gd.T = 4000.0;
        gm.update_trans_coeffs(gd);
        assert(approxEqualNumbers(to!number(0.00012591), gd.mu, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(0.2448263), gd.k, 1.0e-6), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = gd.u; // copy unperturbed value, but we don't really need it
            double ih = 1.0e-20;
            gd.T += complex(0.0,ih);
            gm.update_thermo_from_rhoT(gd);
            double myCv = gd.u.im/ih;
            assert(isClose(myCv, gm.dudT_const_v(gd).re), failedUnitTest());
        }

        // [TODO]
        // entropy, enthalpy and sound speed tests
        return 0;
    } // end main()
}
