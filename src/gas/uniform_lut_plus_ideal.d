/**
 * uniform_lut_plus_ideal.d
 * Composite gas model for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2019-07-09.
 */

module gas.uniform_lut_plus_ideal;

import gas.gas_model;
import gas.gas_state;
import gas.ideal_gas;
import gas.uniform_lut;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;
import nm.complex;
import nm.number;

class UniformLUTPlusIdealGas: GasModel {
public:

    this(lua_State *L) {
        _n_species = 2;
        _n_modes = 0;
        _species_names.length = 2;
        _species_names[0] = "lut";
        _species_names[1] = "ideal";
        // Bring table to TOS
        lua_getglobal(L, "UniformLUTPlusIdealGas");
        // There are just two file-names to be read from the composite-gas-model file.
        string lut_file = getString(L, -1, "lut_file");
        string ideal_file = getString(L, -1, "ideal_file");
        luaL_dofile(L, lut_file.toStringz());
        lut_gas = new UniformLUT(L);
        Q_lut = new GasState(lut_gas);
        luaL_dofile(L, ideal_file.toStringz());
        ideal_gas = new IdealGas(L);
        Q_ideal = new GasState(ideal_gas);
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "UniformLUTPlusIdealGas =(";
        repr ~= "lut_gas=" ~ to!string(lut_gas);
        repr ~= ", ideal_gas=" ~ to!string(ideal_gas);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q)
    {
        if (Q.T <= 0.0 || Q.p <= 0.0) {
            string msg = "Temperature and/or pressure was negative for update_thermo_from_pT."; 
            throw new GasModelException(msg);
        }
        Q_lut.p = Q.massf[0]*Q.p; Q_lut.T = Q.T;
        lut_gas.update_thermo_from_pT(Q_lut);
        Q_ideal.p = Q.massf[1]*Q.p; Q_ideal.T = Q.T;
        ideal_gas.update_thermo_from_pT(Q_ideal);
        Q.rho = Q_lut.rho + Q_ideal.rho;
        Q.u = Q.massf[0]*Q_lut.u + Q.massf[1]*Q_ideal.u;
    }
    override void update_thermo_from_rhou(GasState Q)
    {
        if (Q.u <= 0.0 || Q.rho <= 0.0) {
            string msg = "Internal energy and/or density was negative for update_thermo_from_rhou."; 
            throw new GasModelException(msg);
        }
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        // Need to determine the temperature at which the individual energies sum
        // to the mixture internal energy.
        number T = Q.u/ideal_gas.dudT_const_v(Q_ideal); // First guess is just ideal gas
        Q_lut.T = T; Q_ideal.T = T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        number u_error = (Q.massf[0]*Q_lut.u + Q.massf[1]*Q_ideal.u) - Q.u;
        debug { writefln("u_error=%g", u_error); } 
        // [FIX-ME] finish iterations
        // Partial pressures just sum to the mixture pressure.
        Q.p = Q_lut.p + Q_ideal.p;
        Q.T = T;
    }
    override void update_thermo_from_rhoT(GasState Q)
    {
        if (Q.T <= 0.0 || Q.rho <= 0.0) {
            string msg = "Temperature and/or density was negative for update_thermo_from_rhoT."; 
            throw new GasModelException(msg);
        }
        // Components just add together.
        Q_lut.T = Q.T;
        Q_lut.rho = Q.massf[0] * Q.rho;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        Q_ideal.T = Q.T;
        Q_ideal.rho = Q.massf[1] * Q.rho;
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        Q.p = Q_lut.p + Q_ideal.p;
        Q.u = Q_lut.u + Q_ideal.u;
    }
    override void update_thermo_from_rhop(GasState Q)
    {
        if (Q.p <= 0.0 || Q.rho <= 0.0) {
            string msg = "Pressure and/or density was negative for update_thermo_from_rhop."; 
            throw new GasModelException(msg);
        }
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        // Need to determine the temperature at which the partial pressures sum
        // to the mixture pressure.
        number T = Q.p/ideal_gas.gas_constant(Q_ideal); // First guess is just ideal gas
        Q_lut.T = T; Q_ideal.T = T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        number p_error = (Q_lut.p + Q_ideal.p) - Q.p;
        debug { writefln("p_error=%g", p_error); } 
        // [FIX-ME] finish iterations
        Q.T = T;
        Q.u = Q.massf[0]*Q_lut.u + Q.massf[1]*Q_ideal.u;
    }
    
    override void update_thermo_from_ps(GasState Q, number s)
    {
        // Q.T = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
        // update_thermo_from_pT(Q);
        throw new Error("not implemented");
    }
    override void update_thermo_from_hs(GasState Q, number h, number s)
    {
        // Q.T = h / _Cp;
        // Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.T/_T1)));
        // update_thermo_from_pT(Q);
        throw new Error("not implemented");
    }
    override void update_sound_speed(GasState Q)
    {
        if (Q.T <= 0.0) {
            string msg = "Temperature was negative for update_sound_speed."; 
            throw new GasModelException(msg);
        }
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        lut_gas.update_sound_speed(Q_lut);
        ideal_gas.update_sound_speed(Q_ideal);
        Q.a = Q.massf[0]*Q_lut.a + Q.massf[1]*Q_ideal.a;
    }
    override void update_trans_coeffs(GasState Q)
    {
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        lut_gas.update_trans_coeffs(Q_lut);
        ideal_gas.update_trans_coeffs(Q_ideal);
        Q.mu = Q.massf[0]*Q_lut.mu + Q.massf[1]*Q_ideal.mu;
        Q.k = Q.massf[0]*Q_lut.k + Q.massf[1]*Q_ideal.k;
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
        throw new Exception("not implemented");
    }
    */
    override number dudT_const_v(in GasState Q)
    {
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        return Q.massf[0]*lut_gas.dudT_const_v(Q_lut) +
            Q.massf[1]*ideal_gas.dudT_const_v(Q_ideal);
    }
    override number dhdT_const_p(in GasState Q)
    {
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        return Q.massf[0]*lut_gas.dhdT_const_p(Q_lut) +
            Q.massf[1]*ideal_gas.dhdT_const_p(Q_ideal);
    }
    override number dpdrho_const_T(in GasState Q)
    {
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        return Q.massf[0]*lut_gas.dpdrho_const_T(Q_lut) +
            Q.massf[1]*ideal_gas.dpdrho_const_T(Q_ideal);
    }
    override number gas_constant(in GasState Q)
    {
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        return Q.massf[0]*lut_gas.gas_constant(Q_lut) +
            Q.massf[1]*ideal_gas.gas_constant(Q_ideal);
    }
    override @nogc number internal_energy(in GasState Q)
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q)
    {
        return Q.u + Q.p/Q.rho;
    }
    override number entropy(in GasState Q)
    {
        Q_lut.rho = Q.massf[0]*Q.rho;
        Q_ideal.rho = Q.massf[1]*Q.rho;
        Q_lut.T = Q.T; Q_ideal.T = Q.T;
        lut_gas.update_thermo_from_rhoT(Q_lut);
        ideal_gas.update_thermo_from_rhoT(Q_ideal);
        return Q.massf[0]*lut_gas.entropy(Q_lut) +
            Q.massf[1]*ideal_gas.entropy(Q_ideal);
    }

private:
    // This gas model is composed to two other gas models.
    GasModel lut_gas;
    GasModel ideal_gas;
    GasState Q_lut, Q_ideal;
} // end class UniformLUTPlusIdealGas

version(uniform_lut_plus_ideal_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "uniform-lut-plus-ideal-air-gas-model.lua");
        auto gm = new UniformLUTPlusIdealGas(L);
        lua_close(L);
        auto gd = new GasState(2, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 0.5;
        gd.massf[1] = 0.5;
        assert(approxEqualNumbers(gm.R(gd), to!number(287.086), 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 2, failedUnitTest());
        assert(approxEqualNumbers(gd.p, to!number(1.0e5), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[0], to!number(0.5), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[1], to!number(0.5), 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        writeln("rho=", gd.rho);
        assert(approxEqualNumbers(gd.rho, to!number(1.16113), 1.0e-4), failedUnitTest());
        writeln("u=", gd.u);
        assert(approxEqualNumbers(gd.u, to!number(215643.0), 1.0e-4), failedUnitTest());
        writeln("a=", gd.a);
        assert(approxEqualNumbers(gd.a, to!number(346.625), 1.0e-4), failedUnitTest());
        gm.update_trans_coeffs(gd);
        writeln("mu=", gd.mu);
        assert(approxEqualNumbers(gd.mu, to!number(2.22104e-05), 1.0e-6), failedUnitTest());
        writeln("k=", gd.k);
        assert(approxEqualNumbers(gd.k, to!number(0.0316054), 1.0e-6), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = gd.u; // copy unperturbed value, but we don't really need it
            double h = 1.0e-20;
            gd.T += complex(0.0,h);
            gm.update_thermo_from_rhoT(gd);
            double myCv = gd.u.im/h;
            assert(approxEqual(myCv, gm.dudT_const_v(gd).re), failedUnitTest());
        }
        return 0;
    }
}
