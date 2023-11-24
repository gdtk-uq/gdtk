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
import ntypes.complex;
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
        create_species_reverse_lookup();
        //
        // Underlying thermally-perfect gas model and equilibrium reactor.
        string tpGasEqFile = getString(L, -1, "tpGasEqFile");
        tpGasEqModel = new ThermallyPerfectGasEquilibrium(tpGasEqFile);
        tpgs = GasState(tpGasEqModel);
        tpEqCalc = new EquilibriumCalculator(tpGasEqFile);
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

    GasState savedGasState()
    {
        return tpgs;
    }

    ThermallyPerfectGasEquilibrium savedGasModel()
    {
        return tpGasEqModel;
    }

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

    override void update_thermo_from_pT(ref GasState Q)
    {
        copy_gas_state_except_massf(Q, tpgs);
        tpgs.massf[] = reactants_massf[];

        tpEqCalc.set_massf_from_pT(tpgs);
        tpGasEqModel.update_thermo_from_pT(tpgs);

        copy_gas_state_except_massf(tpgs, Q);
        number Mmix = gas.gas_model.mixture_molecular_mass(tpgs.massf, tpGasEqModel.mol_masses);
        version(complex_numbers) {_mol_masses[0] = Mmix.re;}
        else {_mol_masses[0] = Mmix;}
    }
    override void update_thermo_from_rhou(ref GasState Q)
    {
        copy_gas_state_except_massf(Q, tpgs);
        tpgs.massf[] = reactants_massf[];

        tpEqCalc.set_massf_and_T_from_rhou(tpgs);
        tpGasEqModel.update_thermo_from_rhoT(tpgs); // Faster than calling from_rhou, and ceq set T for us

        copy_gas_state_except_massf(tpgs, Q);
        number Mmix = gas.gas_model.mixture_molecular_mass(tpgs.massf, tpGasEqModel.mol_masses);
        version(complex_numbers) {_mol_masses[0] = Mmix.re;}
        else {_mol_masses[0] = Mmix;}
    }
    override void update_thermo_from_rhoT(ref GasState Q)
    {
        copy_gas_state_except_massf(Q, tpgs);
        tpgs.massf[] = reactants_massf[];

        tpEqCalc.set_massf_from_rhoT(tpgs);
        tpGasEqModel.update_thermo_from_rhoT(tpgs);

        copy_gas_state_except_massf(tpgs, Q);
        number Mmix = gas.gas_model.mixture_molecular_mass(tpgs.massf, tpGasEqModel.mol_masses);
        version(complex_numbers) {_mol_masses[0] = Mmix.re;}
        else {_mol_masses[0] = Mmix;}
    }
    override void update_thermo_from_rhop(ref GasState Q)
    {
        //tpGasEqModel.update_thermo_from_rhop(tpgs); // [TODO]
        throw new Error("Not Implemented Error (equilibrium_gas): update_thermo_from_rhop");
    }
    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        copy_gas_state_except_massf(Q, tpgs);
        tpgs.massf[] = reactants_massf[];

        tpEqCalc.set_massf_and_T_from_ps(tpgs, s.re);
        tpGasEqModel.update_thermo_from_pT(tpgs); // p and T are now set, so call the pT method
        tpGasEqModel.update_trans_coeffs(tpgs);
        tpGasEqModel.update_sound_speed(tpgs);

        copy_gas_state_except_massf(tpgs, Q);
        number Mmix = gas.gas_model.mixture_molecular_mass(tpgs.massf, tpGasEqModel.mol_masses);
        version(complex_numbers) {_mol_masses[0] = Mmix.re;}
        else {_mol_masses[0] = Mmix;}
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        //tpGasEqModel.update_thermo_from_hs(tpgs, h, s);
        throw new Error("Not Implemented Error (equilibrium_gas): update_thermo_from_hs");
    }

    override void update_sound_speed(ref GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        tpGasEqModel.update_sound_speed(tpgs);
        Q.a = tpgs.a;
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        tpGasEqModel.update_trans_coeffs(tpgs);
        Q.mu = tpgs.mu;
        Q.k = tpgs.k;
    }
    override number dudT_const_v(in GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        return tpGasEqModel.dudT_const_v(tpgs);
    }
    override number dhdT_const_p(in GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        return tpGasEqModel.dhdT_const_p(tpgs);
    }
    override number dpdrho_const_T(in GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        return tpGasEqModel.dpdrho_const_T(tpgs);
    }
    override number gas_constant(in GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        return tpGasEqModel.gas_constant(tpgs);
    }
    override number internal_energy(in GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        return tpGasEqModel.internal_energy(tpgs);
    }
    override number enthalpy(in GasState Q)
    {
        set_tpgs_from_external_gasstate(Q);
        return tpGasEqModel.enthalpy(tpgs);
    }
    override number entropy(in GasState Q)
    {
        // It's important that this method return one consistent with ceq's internal thermo routines,
        // so we use those here rather than tpGasEqModel. For some reason tpGasModel has a different
        // zero point to thermo.c, with regard to entropy.
        set_tpgs_from_external_gasstate(Q);
        number s = to!number(tpEqCalc.get_s(tpgs));
        return s;
    }

private:
    ThermallyPerfectGasEquilibrium tpGasEqModel;
    EquilibriumCalculator tpEqCalc;
    number[] reactants_molef;
    number[] reactants_massf;
    GasState tpgs;

    @nogc
    void copy_gas_state_except_massf(ref const(GasState) Q2, ref GasState Q1){
        Q1.rho = Q2.rho;
        Q1.p = Q2.p;
        Q1.p_e = Q2.p_e;
        Q1.T = Q2.T;
        Q1.u = Q2.u;
        Q1.a = Q2.a;
        Q1.mu = Q2.mu;
        Q1.k = Q2.k;
        Q1.sigma = Q2.sigma;
        Q1.quality = Q2.quality;
    }

    @nogc void set_tpgs_from_external_gasstate(ref const(GasState) Q){
        /*
            This is a bandaid to deal with the fact that we have no storage for the
            mass fractions. Instead we use this routine to frantically recalculate them
            whenever the outside world asks this gas model for something.

            Because of this, you shouldn't use this model in performance critical
            calculations (not yet at least...)
        */
        copy_gas_state_except_massf(Q, tpgs);
        tpgs.massf[] = reactants_massf[];
        tpEqCalc.set_massf_from_pT(tpgs);
    }
} // end class EquilibriumGas

// Unit test of the basic gas model...

version(equilibrium_gas_test) {
    import std.stdio;
    import util.msg_service;
    import std.process;

    int main() {
        // Before running this test, we need the gas model files in place.
        auto cmd = executeShell("cp sample-data/air-5sp-eq.lua .");
        assert(cmd.status == 0, failedUnitTest());
        cmd = executeShell("cp ../../examples/kinetics/air-chemistry-1T/air-5sp-1T.inp .");
        assert(cmd.status == 0, failedUnitTest());
        cmd = executeShell("prep-gas air-5sp-1T.inp air-5sp-1T.lua");
        assert(cmd.status == 0, failedUnitTest());

        auto gm = new EquilibriumGas("air-5sp-eq.lua");
        // writeln("gm=", gm); // Can see that the reactants_massf and _molef set.

        auto gd = GasState(1, 0);
        gd.rho = 0.0139284858607; // Fixed example to use reactants = {N2=0.79, O2=0.21} in moles
        gd.u = 2122510.049202302; //

        gd.T = 3352.068185; // Same guess as normal ceq
        gm.update_thermo_from_rhou(gd);
        // writeln("gd.p=", gd.p, " gd.T=", gd.T);
        assert(isClose(2500.0, gd.T, 1.0), failedUnitTest());
        assert(isClose(10135.0, gd.p, 1.0), failedUnitTest());
        return 0;
    }
}
