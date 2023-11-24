/**
 * ideal_dissociating_gas.d
 *
 * This is the Lighthill-Freeman ideal dissociating gas model as described in
 * M. N. Macrossan (1990)
 * Hypervelocity flow of dissociating nitrogen downstream of a blunt nose.
 * Journal of Fluid Mechanics Vol. 217, pp 167-202.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2017-04-22: initial shell copied from powers-aslam-gas module.
 *          2017-04-11: Filled in IDG thermo details as noted in
 *                      PJ's workbooks 2017-04-22 through 2017-05-11
 *          2017-07-27: Finished reactor function and eliminated bugs in workbook.
 */

module gas.ideal_dissociating_gas;

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
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;

// First, the basic gas model.

class IdealDissociatingGas: GasModel {
public:

    this(lua_State *L) {
        type_str = "IdealDissociatingGas";
        // Some parameters are fixed and some come from the gas model file.
        _n_species = 2;
        _n_modes = 0;
        _species_names.length = 2;
        // Bring table to TOS
        lua_getglobal(L, "IdealDissociatingGas");
        // [TODO] test that we actually have the table as item -1
        // Now, pull out the remaining value parameters.
        _species_names[0] = getString(L, -1, "molecule");
        _species_names[1] = getString(L, -1, "atom");
        _W = getDouble(L, -1, "W"); // molecular weight, g/mole
        _Rnn = R_universal / _W * 1000; // gas constant for molecule, J/kg/K
        _mol_masses.length = 2;
        _mol_masses[0] = _W/1000; // kg/mole
        _mol_masses[1] = _W/2000; // kg/mole
        _T_d = getDouble(L, -1, "T_d"); // characteristic dissociation temperature, K
        _rho_d = getDouble(L, -1, "rho_d"); // characteristic density, g/cm^3
        // Rate constants follow.
        _C1 = getDouble(L, -1, "C1");
        _n1 = getDouble(L, -1, "n1");
        _C2 = getDouble(L, -1, "C2");
        _n2 = getDouble(L, -1, "n2");
        lua_pop(L, 1); // dispose of the table
        create_species_reverse_lookup();
        // Reference conditions for entropy.
        _Tref = 298.15; _pref = 101.325e3;
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= "IdealDissociatingGas =(";
        repr ~= "species=[\""~_species_names[0];
        repr ~= "\", \""~_species_names[1]~"\"]";
        repr ~= ", W=", to!string(_W);
        repr ~= ", Mmass=[" ~ to!string(_mol_masses[0]);
        repr ~= "," ~ to!string(_mol_masses[1]) ~ "]";
        repr ~= ", T_d=" ~ to!string(_T_d);
        repr ~= ", rho_d=" ~ to!string(_rho_d);
        repr ~= ", C1=" ~ to!string(_C1);
        repr ~= ", n1=" ~ to!string(_n1);
        repr ~= ", C2=" ~ to!string(_C2);
        repr ~= ", n2=" ~ to!string(_n2);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q) const
    {
        number alpha = Q.massf[1];
        Q.rho = Q.p/(Q.T*_Rnn*(1+alpha));
        Q.u = _Rnn*alpha*_T_d + _Rnn*3*Q.T;
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        number alpha = Q.massf[1];
        Q.T = (Q.u - _Rnn*alpha*_T_d)/(_Rnn*3);
        Q.p = Q.rho*(1+alpha)*_Rnn*Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        number alpha = Q.massf[1];
        Q.p = Q.rho*(1+alpha)*_Rnn*Q.T;
        Q.u = _Rnn*alpha*_T_d + _Rnn*3*Q.T;
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        number alpha = Q.massf[1];
        Q.T = Q.p/(Q.rho*(1+alpha)*_Rnn*Q.T);
        Q.u = _Rnn*alpha*_T_d + _Rnn*3*Q.T;
    }

    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        // For frozen composition.
        number alpha = Q.massf[1];
        Q.T = _Tref * exp((1.0/((4+alpha)*_Rnn))*(s+(1+alpha)*_Rnn*log(Q.p/_pref)));
        update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        // For frozen composition.
        number alpha = Q.massf[1];
        Q.T = (h - _Rnn*alpha*_T_d) / ((4+alpha)*_Rnn);
        Q.p = _pref * exp(((4+alpha)*_Rnn*log(Q.T/_Tref)-s)/((4+alpha)*_Rnn));
        update_thermo_from_pT(Q);
    }
    override void update_sound_speed(ref GasState Q) const
    {
        // For frozen composition.
        number alpha = Q.massf[1];
        number gamma = (4+alpha)/3.0;
        number Rgas = (1+alpha)*_Rnn;
        Q.a = sqrt(gamma*Rgas*Q.T);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        // The gas is inviscid.
        Q.mu = 0.0;
        Q.k = 0.0;
    }
    override number dudT_const_v(in GasState Q) const
    {
        number alpha = Q.massf[1];
        return to!number(_Rnn*3); // frozen alpha
    }
    override number dhdT_const_p(in GasState Q) const
    {
        number alpha = Q.massf[1];
        return _Rnn*(4+alpha); // frozen alpha
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        number alpha = Q.massf[1];
        return _Rnn*(1+alpha)*Q.T; // frozen alpha
    }
    override number gas_constant(in GasState Q) const
    {
        number alpha = Q.massf[1];
        return _Rnn*(1+alpha);
    }
    override number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q) const
    {
        return Q.u + Q.p/Q.rho;
    }
    override number entropy(in GasState Q) const
    {
        // Presume that we have a fixed composition mixture and that
        // Entropy for each species is zero at reference condition.
        number alpha = Q.massf[1];
        return (4+alpha)*log(Q.T/_Tref) - (1+alpha)*_Rnn*log(Q.p/_pref);
    }

private:
    // Thermodynamic constants
    double _W; // Molecular mass in g/mole
    double _Rnn; // gas constant for molecule in J/kg/K
    double _T_d; // characteristic dissociation temperature, K
    double _rho_d; // characteristic density, g/cm^3
    double _Tref, _pref; // Reference conditions for entropy
    // Rate constants
    double _C1, _n1, _C2, _n2;
   // Molecular transport coefficents are zero.
} // end class IdealDissociatingGas


// Unit test of the basic gas model...

version(ideal_dissociating_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/idg-nitrogen.lua");
        auto gm = new IdealDissociatingGas(L);
        lua_close(L);
        auto gd = GasState(2, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 1.0; gd.massf[1] = 0.0;
        number Rgas = R_universal / 0.028;
        assert(isClose(gm.R(gd), Rgas, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 2, failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[1], 0.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        number my_rho = 1.0e5 / (Rgas * 300.0);
        assert(isClose(gd.rho, my_rho, 1.0e-4), failedUnitTest());
        number my_Cv = gm.dudT_const_v(gd);
        number my_u = my_Cv*300.0;
        assert(isClose(gd.u, my_u, 1.0e-3), failedUnitTest());
        number my_Cp = gm.dhdT_const_p(gd);
        number my_a = sqrt(my_Cp/my_Cv*Rgas*300.0);
        assert(isClose(gd.a, my_a, 1.0e-3), failedUnitTest());
        gm.update_trans_coeffs(gd);
        assert(isClose(gd.mu, 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.k, 0.0, 1.0e-6), failedUnitTest());

        // Let's set a higher temperature condition so that we can see some reaction.
        // [TODO]: move this test to kinetics/ideal_dissociating_gas_kinetics.d
        /*
        gd.T = 5000.0;
        gd.p = 1.0e5;
        gm.update_thermo_from_pT(gd);
        auto reactor = new UpdateIDG("sample-data/idg-nitrogen.lua", gm);
        number[] params; // empty
        number dtSuggest; // to receive value
        reactor(gd, 1.0e-3, dtSuggest, params);
        */

        return 0;
    } // end main
}
