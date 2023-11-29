/**
 * cubic_gas.d
 * Dense-gas model with a cubic EOS for use in the CFD codes.
 *
 * Author: Matthew Kratzer, Peter J. and Rowan G.
 * Version: 2019-05-10: initial cut adapted from ideal_gas.
 */

module gas.cubic_gas;

import gas.gas_model;
import gas.gas_state;
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
import ntypes.complex;
import nm.number;

class CubicGas: GasModel {
public:

    this(lua_State *L) {
        type_str = "CubicGas";
        _n_species = 1;
        _n_modes = 0;
        // Bring table to TOS
        lua_getglobal(L, "CubicGas");
        // There's just one species name to be read from the gas-model file.
        _species_names.length = 1;
        _species_names[0] = getString(L, -1, "speciesName");
        // Now, pull out the remaining numeric value parameters.
        _mol_masses.length = 1;
        _mol_masses[0] = getDouble(L, -1, "mMass");
        _gamma = getDouble(L, -1, "gamma");
        // Reference values for entropy
        lua_getfield(L, -1, "entropyRefValues");
        _s1 = getDouble(L, -1, "s1");
        _T1 = getDouble(L, -1, "T1");
        _p1 = getDouble(L, -1, "p1");
        lua_pop(L, 1);
        // Molecular transport coefficent models.
        lua_getfield(L, -1, "viscosity");
        _viscModel = createViscosityModel(L);
        lua_pop(L, 1);

        lua_getfield(L, -1, "thermCondModel");
        auto model = getString(L, -1, "model");
        if ( model == "constPrandtl" ) {
            _Prandtl = getDouble(L, -1, "Prandtl");
            _constPrandtl = true;
        }
        else {
            _thermCondModel = createThermalConductivityModel(L);
            _constPrandtl = false;
        }
        lua_pop(L, 1);
        // Compute derived parameters
        _Rgas = R_universal/_mol_masses[0];
        _Cv = _Rgas / (_gamma - 1.0);
        _Cp = _Rgas*_gamma/(_gamma - 1.0);
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "CubicGas =(";
        repr ~= "name=\"" ~ _species_names[0] ~"\"";
        repr ~= ", Mmass=" ~ to!string(_mol_masses[0]);
        repr ~= ", gamma=" ~ to!string(_gamma);
        repr ~= ", s1=" ~ to!string(_s1);
        repr ~= ", T1=" ~ to!string(_T1);
        repr ~= ", p1=" ~ to!string(_p1);
        repr ~= ", constPrandtl=" ~ to!string(_constPrandtl);
        repr ~= ", Prandtl=" ~ to!string(_Prandtl);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q) const
    {
        if (Q.T <= 0.0 || Q.p <= 0.0) {
            string msg = "Temperature and/or pressure was negative for update_thermo_from_pT.";
            throw new GasModelException(msg);
        }
        Q.rho = Q.p/(Q.T*_Rgas);
        Q.u = _Cv*Q.T;
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        if (Q.u <= 0.0 || Q.rho <= 0.0) {
            string msg = "Internal energy and/or density was negative for update_thermo_from_rhou.";
            throw new GasModelException(msg);
        }
        Q.T = Q.u/_Cv;
        Q.p = Q.rho*_Rgas*Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        if (Q.T <= 0.0 || Q.rho <= 0.0) {
            string msg = "Temperature and/or density was negative for update_thermo_from_rhoT.";
            throw new GasModelException(msg);
        }
        Q.p = Q.rho*_Rgas*Q.T;
        Q.u = _Cv*Q.T;
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        if (Q.p <= 0.0 || Q.rho <= 0.0) {
            string msg = "Pressure and/or density was negative for update_thermo_from_rhop.";
            throw new GasModelException(msg);
        }
        Q.T = Q.p/(Q.rho*_Rgas);
        Q.u = _Cv*Q.T;
    }

    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        Q.T = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
        update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        Q.T = h / _Cp;
        Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.T/_T1)));
        update_thermo_from_pT(Q);
    }
    override void update_sound_speed(ref GasState Q) const
    {
        if (Q.T <= 0.0) {
            string msg = "Temperature was negative for update_sound_speed.";
            throw new GasModelException(msg);
        }
        Q.a = sqrt(_gamma*_Rgas*Q.T);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        Q.mu = _viscModel.eval(Q.T);
        if ( _constPrandtl ) {
            Q.k = _Cp*Q.mu/_Prandtl;
        }
        else {
            Q.k = _thermCondModel.eval(Q.T);
        }
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
        throw new Exception("not implemented");
    }
    */
    override number dudT_const_v(in GasState Q) const
    {
        return to!number(_Cv);
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return to!number(_Cp);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        number R = gas_constant(Q);
        return R*Q.T;
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(R_universal/_mol_masses[0]);
    }
    override @nogc number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q) const
    {
        return Q.u + Q.p/Q.rho;
    }
    override number entropy(in GasState Q) const
    {
        return _s1 + _Cp * log(Q.T/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    // Thermodynamic constants
    double _Rgas; // J/kg/K
    double _gamma;   // ratio of specific heats
    double _Cv; // J/kg/K
    double _Cp; // J/kg/K
    // Reference values for entropy
    double _s1;  // J/kg/K
    double _T1;  // K
    double _p1;  // Pa
    // Molecular transport coefficent constants.
    Viscosity _viscModel;
    // We compute thermal conductivity in one of two ways:
    // 1. based on constant Prandtl number; OR
    // 2. ThermalConductivity model
    // therefore we have places for both data
    bool _constPrandtl = false;
    double _Prandtl;
    ThermalConductivity _thermCondModel;

} // end class CubicGas

version(cubic_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/cubic-gas-model.lua");
        auto gm = new CubicGas(L);
        lua_close(L);
        auto gd = GasState(1, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 1.0;
        assert(approxEqualNumbers(gm.R(gd), to!number(287.086), 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(approxEqualNumbers(gd.p, to!number(1.0e5), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.massf[0], to!number(1.0), 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        assert(approxEqualNumbers(gd.rho, to!number(1.16109), 1.0e-4), failedUnitTest());
        assert(approxEqualNumbers(gd.u, to!number(215314.0), 1.0e-4), failedUnitTest());
        assert(approxEqualNumbers(gd.a, to!number(347.241), 1.0e-4), failedUnitTest());
        gm.update_trans_coeffs(gd);
        assert(approxEqualNumbers(gd.mu, to!number(1.84691e-05), 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(gd.k, to!number(0.0262449), 1.0e-6), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = gd.u; // copy unperturbed value, but we don't really need it
            double h = 1.0e-20;
            gd.T += complex(0.0,h);
            gm.update_thermo_from_rhoT(gd);
            double myCv = gd.u.im/h;
            assert(isClose(myCv, gm.dudT_const_v(gd).re), failedUnitTest());
        }
        return 0;
    }
}
