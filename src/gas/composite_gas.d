/**
 * composite_gas.d
 * A gas model with behaviour built by composition of pieces.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2021-02-09
 */

module gas.composite_gas;

import std.string;
import std.conv;
import std.algorithm;

import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;

import gas.thermo.therm_perf_gas_mix;
import gas.diffusion.gas_mixtures;

class CompositeGas : GasModel {
public:

    this(lua_State *L)
    {
        type_str = "CompositeGas";
        /* There are some top-level GasModel services that require us to fill in data
         * at this point.
         */
        getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", _species_names);
        _n_species = cast(uint) _species_names.length;
        create_species_reverse_lookup();
        _n_modes = 0; // Single temperature gas
        if (canFind(_species_names, "e-") || canFind(_species_names, "eminus")) {
            if (!(_species_names[$-1] == "e-" || _species_names[$-1] == "eminus")) {
                throw new Error("Electrons should be last species.");
            }
            _is_plasma = true;
        }

        _mol_masses.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _mol_masses[isp] = getDouble(L, -1, "M");
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        mThermo = new ThermPerfGasMixture(L, _species_names);
        mTransProps = new GasMixtureTransProps(L, _species_names);

        // Fill in some parameters needed in GasModel protected data.
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
    }
    
    // Service methods related to thermodynamics
    // Updates to GasState
    @nogc override void update_thermo_from_pT(GasState gs)
    {
        mThermo.updateFromPT(gs);
    }
    @nogc override void update_thermo_from_rhou(GasState gs)
    {
        mThermo.updateFromRhoU(gs);
    }
    @nogc override void update_thermo_from_rhoT(GasState gs)
    {
        mThermo.updateFromRhoT(gs);
    }
    @nogc override void update_thermo_from_rhop(GasState gs)
    {
        mThermo.updateFromRhoP(gs);
    }
    @nogc override void update_thermo_from_ps(GasState gs, number s)
    {
        mThermo.updateFromPS(gs, s);
    }
    @nogc override void update_thermo_from_hs(GasState gs, number h, number s)
    {
        mThermo.updateFromHS(gs, h, s);
    }
    @nogc override void update_sound_speed(GasState gs)
    {
        mThermo.updateSoundSpeed(gs);
    }
    // Return of single values
    @nogc override number dudT_const_v(in GasState gs)
    {
        return mThermo.dudTConstV(gs);
    }
    @nogc override number dhdT_const_p(in GasState gs)
    {
        return mThermo.dhdTConstP(gs);
    }
    @nogc override number dpdrho_const_T(in GasState gs)
    {
        return mThermo.dpdrhoConstT(gs);
    }
    @nogc override number gas_constant(in GasState gs)
    {
        return mThermo.gasConstant(gs);
    }
    @nogc override number internal_energy(in GasState gs)
    {
        return mThermo.internalEnergy(gs);
    }
    @nogc override number enthalpy(in GasState gs)
    {
        return mThermo.enthalpy(gs);
    }
    @nogc override number enthalpy(in GasState gs, int isp)
    {
        return mThermo.enthalpy(gs, isp);
    }
    @nogc override number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        return mThermo.enthalpy(gs, isp, imode);
    }
    @nogc override number entropy(in GasState gs)
    {
        return mThermo.entropy(gs);
    }
    @nogc override number entropy(in GasState gs, int isp)
    {
        return mThermo.entropy(gs, isp);
    }

    @nogc override void update_trans_coeffs(GasState gs)
    {
        mTransProps.updateTransProps(gs);
    }

private:
    ThermPerfGasMixture mThermo;
    GasMixtureTransProps mTransProps;
}


version(composite_gas_test) {
    int main() {
        import util.msg_service;
        import std.math;
        
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
        //
        auto L = init_lua_State();
        doLuaFile(L, "sample-data/therm-perf-5-species-air.lua");
        auto gm = new CompositeGas(L);
        lua_close(L);
        auto gs = new GasState(5, 0);
        assert(approxEqual(3.621, gm.LJ_sigmas[0]), failedUnitTest());
        assert(approxEqual(97.530, gm.LJ_epsilons[0]), failedUnitTest());

        gs.p = 1.0e6;
        gs.T = 2000.0;
        gs.massf = [to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2)];
        gm.update_thermo_from_pT(gs);
        assert(approxEqualNumbers(to!number(11801825.6), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.2840117), gs.rho, 1.0e-6), failedUnitTest());

        gs.rho = 2.0;
        gs.u = 14.0e6;
        gm.update_thermo_from_rhou(gs);
        assert(approxEqualNumbers(to!number(3373757.4), gs.p, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(4331.944), gs.T, 1.0e-6), failedUnitTest());

        gs.T = 10000.0;
        gs.rho = 1.5;
        gm.update_thermo_from_rhoT(gs);
        assert(approxEqualNumbers(to!number(5841068.3), gs.p, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(20340105.9), gs.u, 1.0e-6), failedUnitTest());

        gs.rho = 10.0;
        gs.p = 5.0e6;
        gm.update_thermo_from_rhop(gs);
        assert(approxEqualNumbers(to!number(11164648.5), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1284.012), gs.T, 1.0e-6), failedUnitTest());

        gs.p = 1.0e6;
        number s = 10000.0;
        gm.update_thermo_from_ps(gs, s);
        assert(approxEqualNumbers(to!number(2560.118), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(12313952.52), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.00309), gs.rho, 1.0e-6), failedUnitTest());

        s = 11000.0;
        number h = 17.0e6;
        gm.update_thermo_from_hs(gs, h, s);
        assert(approxEqualNumbers(to!number(5273.103), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(14946629.7), gs.u, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(0.4603513), gs.rho, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(945271.84), gs.p, 1.0e-4), failedUnitTest());

        gs.T = 4000.0;
        gm.update_trans_coeffs(gs);
        assert(approxEqualNumbers(to!number(0.00012591), gs.mu, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(0.2448263), gs.k, 1.0e-6), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = gs.u; // copy unperturbed value, but we don't really need it
            double ih = 1.0e-20;
            gs.T += complex(0.0,ih);
            gm.update_thermo_from_rhoT(gs);
            double myCv = gs.u.im/ih;
            assert(approxEqual(myCv, gm.dudT_const_v(gs).re), failedUnitTest());
        }

        return 0;
    } // end main()

}
