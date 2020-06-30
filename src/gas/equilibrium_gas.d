/**
 * equilibrium_gas.d
 *
 * A thermally-perfect gas mix with equilibrium thermochemistry.
 *
 * Authors: Nick Gibons, Peter J. and Rowan G.
 * Version: 2020-06-30
 */

module gas.equilibrium_gas;

import std.math;
import std.stdio;
import std.string;
import std.algorithm;
import std.conv;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.therm_perf_gas_equil;

import kinetics.thermochemical_reactor;
import kinetics.equilibrium_update;

// First, the basic gas model.

class EquilibriumGas: GasModel {
public:
    this(lua_State *L)
    // Construct model from parameters contained in a Lua interpreter.
    {
        type_str = "EquilibriumGas";
        // In the context of Rowan's perfect gas mix, the equilibrium gas is a strange beast.
        // We will hide it's internal species behind a single pseudo-species that
        // we will call by the mixtureName.  It will not have a fixed molecular mass.
        // Construct from information in a Lua table.
        lua_getglobal(L, "EquilibriumGas"); // Bring that table to TOS
        if (!lua_istable(L, -1)) { throw new Exception("Could not find table named EquilibriumGas."); }
        string mixtureName = getString(L, -1, "mixtureName");
        _n_modes = 0;
        _n_species = 1;
        _species_names.length = _n_species;
        _species_names[0] = mixtureName;
        _mol_masses.length = 1;
        _mol_masses[0] = 0.0; // dummy value; we shouldn't use it
        //
        // Underlying thermally-perfect gas model and equilibrium reactor.
        string tpGasEqFile = getString(L, -1, "tpGasEqFile");
        tpGasEqModel = new ThermallyPerfectGasEquilibrium(tpGasEqFile);
        tpgs = new GasState(tpGasEqModel);
        // [TODO] get rid of reactor and use Nick's ceq functions directly.
        tpGasReactor = new EquilibriumUpdate(tpGasEqFile, tpGasEqModel);
        //
        // For the amounts of the reactants, we are expecting
        // to find only the non-zero components in the Lua file
        // but we will use whatever we find.
        // We will get the species list from the underlying thermally-perfect gas model.
        lua_getfield(L, -1, "reactants");
        double[string] reactants;
        foreach (sname; tpGasEqModel._species_names) {
            lua_getfield(L, -1, sname.toStringz);
            if (lua_isnumber(L, -1)) {
                reactants[sname] = lua_tonumber(L, -1);
            } else {
                reactants[sname] = 0.0;
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1); // dispose of reactants table
        string inputUnits = getString(L, -1, "inputUnits");
        reactants_molef.length = tpGasEqModel.n_species;
        reactants_massf.length = tpGasEqModel.n_species;
        if (inputUnits.canFind("moles")) {
            foreach (i, species; tpGasEqModel._species_names) { reactants_molef[i] = reactants[species]; }
            gas.gas_model.molef2massf(reactants_molef, tpGasEqModel._mol_masses, reactants_massf);
        } else {
            foreach (i, species; tpGasEqModel._species_names) { reactants_massf[i] = reactants[species]; }
            gas.gas_model.massf2molef(reactants_massf, tpGasEqModel._mol_masses, reactants_molef);
        }
    } // end constructor

    this(in string fname)
    // Construct model from a Lua script file.
    {
        auto L = init_lua_State();
        doLuaFile(L, fname);
        this(L);
        lua_close(L);
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= "EquilibriumGas =(";
        repr ~= "reactants_massf=" ~ to!string(reactants_massf);
        repr ~= ", reactants_molef=" ~ to!string(reactants_molef);
        repr ~= ", tpGasEqModel=" ~ to!string(tpGasEqModel);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q)
    {
        tpgs.p = Q.p;
        tpgs.T = Q.T;
        tpgs.massf[] = reactants_massf[];
        // [FIX-ME] this one is not ready to try.
        //
        // [TODO] would be better to use Nick's ceq functions directly,
        // rather than via the reactor.
        // tpGasEqModel.update_thermo_from_pT(tpgs);
        double tInterval = 0.0;
        double dtChemSuggest = 0.0;
        double dtThermSuggest = 0.0;
        number[maxParams] params;
        tpGasReactor(tpgs, tInterval, dtChemSuggest, dtThermSuggest, params);
        Q.rho = tpgs.rho;
        Q.u = tpgs.u;
    }
    override void update_thermo_from_rhou(GasState Q)
    {
        tpgs.rho = Q.rho;
        tpgs.u = Q.u;
        tpgs.massf[] = reactants_massf[];
        // [TODO] would be better to use Nick's ceq functions directly,
        // rather than via the reactor.
        tpgs.T = 2500.0; // To avoid a nan in cea_thermo_curves.d, but should not need it. [FIX-ME]
        tpGasEqModel.update_thermo_from_rhou(tpgs);
        double tInterval = 0.0;
        double dtChemSuggest = 0.0;
        double dtThermSuggest = 0.0;
        number[maxParams] params;
        tpGasReactor(tpgs, tInterval, dtChemSuggest, dtThermSuggest, params);
        Q.p = tpgs.p;
        Q.T = tpgs.T;
    }
    override void update_thermo_from_rhoT(GasState Q)
    {
        tpGasEqModel.update_thermo_from_rhoT(tpgs); // [TODO]
    }
    override void update_thermo_from_rhop(GasState Q)
    {
        tpGasEqModel.update_thermo_from_rhop(tpgs); // [TODO]
    }
    override void update_thermo_from_ps(GasState Q, number s)
    {
        tpGasEqModel.update_thermo_from_ps(tpgs, s);// [TODO]
    }
    override void update_thermo_from_hs(GasState Q, number h, number s)
    {
        tpGasEqModel.update_thermo_from_hs(tpgs, h, s);// [TODO]
    }
    override void update_sound_speed(GasState Q)
    {
        tpGasEqModel.update_sound_speed(tpgs);// [TODO]
    }
    override void update_trans_coeffs(GasState Q)
    {
        tpGasEqModel.update_trans_coeffs(tpgs);// [TODO]
    }
    override number dudT_const_v(in GasState Q)
    {
        return tpGasEqModel.dudT_const_v(tpgs);// [TODO]
    }
    override number dhdT_const_p(in GasState Q)
    {
        return tpGasEqModel.dhdT_const_p(tpgs);// [TODO]
    }
    override number dpdrho_const_T(in GasState Q)
    {
        return tpGasEqModel.dpdrho_const_T(tpgs);// [TODO]
    }
    override number gas_constant(in GasState Q)
    {
        return tpGasEqModel.gas_constant(tpgs);// [TODO]
    }
    override number internal_energy(in GasState Q)
    {
        return tpGasEqModel.internal_energy(tpgs);// [TODO]
    }
    override number enthalpy(in GasState Q)
    {
        return tpGasEqModel.enthalpy(tpgs);// [TODO]
    }
    override number entropy(in GasState Q)
    {
        return tpGasEqModel.entropy(tpgs);// [TODO]
    }

    ThermallyPerfectGasEquilibrium tpGasEqModel;
    EquilibriumUpdate tpGasReactor;
    number[] reactants_molef;
    number[] reactants_massf;
    GasState tpgs;
} // end class EquilibriumGas

// Unit test of the basic gas model...

version(equilibrium_gas_test) {
    import std.stdio;
    import util.msg_service;
    import std.process;

    int main() {
        // Before running this test, we need the gas model files in place.
        auto cmd = executeShell("cp sample-data/air-eq-1.lua .");
        assert(cmd.status == 0, failedUnitTest());
        cmd = executeShell("cp sample-data/air-eq.inp .");
        assert(cmd.status == 0, failedUnitTest());
        cmd = executeShell("prep-gas air-eq.inp air-eq.lua");
        assert(cmd.status == 0, failedUnitTest());

        auto gm = new EquilibriumGas("air-eq-1.lua");
        // writeln("gm=", gm); // Can see that the reactants_massf and _molef set.

        auto gd = new GasState(1, 0);
        gd.rho = 0.0139638507337; // These numbers pulled from Nick's reaction example.
        gd.u = 2131154.032;
        gd.T = 2500.0;
        gm.update_thermo_from_rhou(gd);
        // writeln("gd.p=", gd.p, " gd.T=", gd.T);
        assert(approxEqual(2505.0, gd.T, 1.0), failedUnitTest());
        assert(approxEqual(10182.0, gd.p, 1.0), failedUnitTest());
        return 0;
    }
}
