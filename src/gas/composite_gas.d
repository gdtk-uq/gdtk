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
import std.stdio;

import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;

import gas.thermo.thermo_model;
import gas.thermo.therm_perf_gas_mix;
import gas.thermo.two_temperature_gas;
import gas.thermo.three_temperature_gas;
import gas.thermo.multi_temperature_gas;
import gas.diffusion.transport_properties_model;
import gas.diffusion.gas_mixtures;
import gas.diffusion.two_temperature_trans_props;
import gas.diffusion.three_temperature_trans_props;
import gas.diffusion.multi_temperature_trans_props;


class CompositeGas : GasModel {
public:
    @property string physicalModel() { return mPhysicalModel; }
    this(lua_State *L)
    {
        type_str = "CompositeGas";
        /* There are some top-level GasModel services that require us to fill in data
         * at this point.
         */
        getArrayOfStrings(L, "species", _species_names);
        _n_species = cast(uint) _species_names.length;
        create_species_reverse_lookup();
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

        mPhysicalModel = getString(L, "physical_model");
        switch (mPhysicalModel) {
        case "thermally-perfect-gas":
            _n_modes = 0;
            mThermo = new ThermPerfGasMixture(L, _species_names);
            mTransProps = new GasMixtureTransProps(L, _species_names);
            break;
        case "two-temperature-gas":
            _n_modes = 1;
            _energy_mode_names.length = _n_modes;
            _energy_mode_names[0] = "vibroelectronic";
            create_energy_mode_reverse_lookup();
            mThermo = new TwoTemperatureGasMixture(L, _species_names);
            mTransProps = new TwoTemperatureTransProps(L, _species_names);
            break;
        case "three-temperature-gas":
            _n_modes = 2;
            _energy_mode_names.length = _n_modes;
            _energy_mode_names[0] = "vibrational";
            _energy_mode_names[1] = "electronic";
            create_energy_mode_reverse_lookup();
            mThermo = new ThreeTemperatureGasMixture(L, _species_names);
            mTransProps = new ThreeTemperatureTransProps(L, _species_names);
            break;
        case "multi-temperature-gas":
            getArrayOfStrings(L, "energy_modes", _energy_mode_names);
            _n_modes = to!int(_energy_mode_names.length);
            create_energy_mode_reverse_lookup();
            mThermo = new MultiTemperatureGasMixture(L, _species_names, _energy_mode_names);
            mTransProps = new MultiTemperatureTransProps(L, _species_names, _energy_mode_names);
            break;
        default:
            string errMsg = format("Problem trying to construct gas model. The physical model variant '%s' is not known.\n", mPhysicalModel);
            throw new Error(errMsg);
        }

        // Fill in some parameters needed in GasModel protected data.
        _LJ_sigmas.length = _n_species;
        _LJ_epsilons.length = _n_species;
        _Lewis_numbers.length = _n_species;
        _charge.length = n_species;
        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _LJ_sigmas[isp] = getDouble(L, -1, "sigma");
            _LJ_epsilons[isp] = getDouble(L, -1, "epsilon");
            _Lewis_numbers[isp] = getDouble(L, -1, "Lewis");
            _charge[isp] = getInt(L, -1, "charge");
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
    }

    this(string gas_file_name){
        lua_State* L = init_lua_State();
        doLuaFile(L, gas_file_name);
        this(L);
        lua_close(L);
    }

    // Service methods related to thermodynamics
    // Updates to GasState
    override void update_thermo_from_pu(ref GasState gs)
    {
	mThermo.updateFromPU(gs);
    }
    @nogc override void update_thermo_from_pT(ref GasState gs)
    {
        mThermo.updateFromPT(gs);
    }
    @nogc override void update_thermo_from_rhou(ref GasState gs)
    {
        mThermo.updateFromRhoU(gs);
    }
    @nogc override void update_thermo_from_rhoT(ref GasState gs)
    {
        mThermo.updateFromRhoT(gs);
    }
    @nogc override void update_thermo_from_rhop(ref GasState gs)
    {
        mThermo.updateFromRhoP(gs);
    }
    @nogc override void update_thermo_from_ps(ref GasState gs, number s)
    {
        mThermo.updateFromPS(gs, s);
    }
    @nogc override void update_thermo_from_hs(ref GasState gs, number h, number s)
    {
        mThermo.updateFromHS(gs, h, s);
    }
    @nogc override void update_sound_speed(ref GasState gs)
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
    @nogc override number gas_constant(in GasState gs, int isp)
    {
        return mThermo.gasConstant(gs, isp);
    }
    @nogc override number internal_energy(in GasState gs)
    {
        return mThermo.internalEnergy(gs);
    }
    @nogc override number energyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        return mThermo.energyPerSpeciesInMode(gs, isp, imode);
    }
    @nogc override number enthalpy(in GasState gs)
    {
        return mThermo.enthalpy(gs);
    }
    @nogc override void enthalpies(in GasState gs, number[] hs)
    {
        return mThermo.enthalpies(gs, hs);
    }
    @nogc override number enthalpy(in GasState gs, int isp)
    {
        return mThermo.enthalpyPerSpecies(gs, isp);
    }
    @nogc override number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        return mThermo.enthalpyPerSpeciesInMode(gs, isp, imode);
    }
    @nogc override number entropy(in GasState gs)
    {
        return mThermo.entropy(gs);
    }
    @nogc override number entropy(in GasState gs, int isp)
    {
        return mThermo.entropyPerSpecies(gs, isp);
    }
    @nogc override void gibbs_free_energies(ref GasState Q, number[] gibbs_energies)
    {
        mThermo.GibbsFreeEnergies(Q, gibbs_energies);
    }
    @nogc override number Cp(in GasState Q, int isp)
    {
        return mThermo.cpPerSpecies(Q, isp);
    }

    @nogc override void update_trans_coeffs(ref GasState gs)
    {
        mTransProps.updateTransProps(gs, this);
    }

    @nogc override void binary_diffusion_coefficients(ref const(GasState) Q, ref number[][] D)
    {
        debug{ assert(D.length==_n_species); }
        mTransProps.binaryDiffusionCoefficients(Q, D);
    }

    @nogc override void minimum_mixture_energy(ref GasState Q)
    {
        Q.p = P_atm;
        Q.T = T_MIN;
        foreach (ref T; Q.T_modes) T = T_MIN;
        auto uTot = mThermo.internalEnergy(Q);
        foreach (imode; 0 .. _n_modes) {
            Q.u_modes[imode] = 0.0;
            foreach (isp; 0 .. _n_species) {
                Q.u_modes[imode] += Q.massf[isp] * mThermo.energyPerSpeciesInMode(Q, isp, imode);
            }
        }
        Q.u = uTot - sum(Q.u_modes);
    }

private:
    string mPhysicalModel;
    ThermodynamicModel mThermo;
    TransportPropertiesModel mTransProps;
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
        auto gs = GasState(5, 0);
        assert(isClose(3.621, gm.LJ_sigmas[0]), failedUnitTest());
        assert(isClose(97.530, gm.LJ_epsilons[0]), failedUnitTest());

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
            assert(isClose(myCv, gm.dudT_const_v(gs).re), failedUnitTest());
        }

        return 0;
    } // end main()

}
